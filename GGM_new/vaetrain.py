import argparse
from collections import OrderedDict
import os
import random
import re
import time

import numpy as np
import torch
import torch.multiprocessing as mp
from torch.multiprocessing import Pool
import torch.nn as nn
import torch.optim as optim

from ggm import ggm
from shared_optim import SharedRMSprop, SharedAdam
import utils

def train(shared_model, optimizer, wholes, scaffolds, whole_conditions, scaffold_conditions, pid, retval_list, args):
    """\
    Target function for the multiprocessed training.

    In addition to updating model parameters, 
    loss values are collected by `retval_list` after each `forward`.

    Parameters
    ----------
    shared_model: torch.nn.Module
        A shared model to be trained.
    optimizer: torch.optim.Optimizer
        A shared optimizer.
    wholes: list[str]
        A list of whole-molecule SMILESs.
    scaffolds: list[str]
        A list of scaffold SMILESs.
    whole_conditions: list[ list[float] ]
        [
            [ value1, value2, ... ],  # condition values of whole 1
            [ value1, value2, ... ],  # condition values of whole 2
        ]
    scaffold_conditions: list[ list[float] ]
        Similar to `whole_conditions`, but with scaffold values.
    pid: int
        CPU index.
    retval_list: list[multiprocessing.managers.ListProxy]
        A list of lists to collect loss floats from CPUs.
        In each cycle, the final shape will be:
            (ncpus, minibatch_size, num_of_losses)
    args: argparse.Namespace
        Delivers parameters from command arguments to the model.
    """
    #each thread make new model


    model=ggm(args)


    for idx in range(len(wholes)):
        #set parameters of model as same as that of reference model



        model.load_state_dict(shared_model.state_dict())
        model.zero_grad()
        optimizer.zero_grad()
        
        #forward
        retval = model(wholes[idx], scaffolds[idx], whole_conditions[idx], scaffold_conditions[idx], args.shuffle_order)
        
        #if retval is None, some error occured. it is usually due to invalid smiles
        if retval is None:
            print('RETVAL IS NONE....PANIC!')
            continue

        #train model
        g_gen, h_gen, loss1, loss2, loss3 = retval
        
        #print(f'loss1: {loss1}, {type(loss1)}')
        #print(f'loss2: {loss2}, {type(loss2)}')
        #print(f'loss3: {loss3}, {type(loss3)}')


        loss = loss1 + loss2*args.beta1 + loss3  # torch.autograd.Variable of shape (1,)
        #print(loss.data.cpu().numpy())

        #retval_list[pid].append((loss.data.cpu().numpy()[0], loss1.data.cpu().numpy()[0], loss2.data.cpu().numpy()[0], loss3.data.cpu().numpy()[0]))
        retval_list[pid].append((loss.data.cpu().numpy(), loss1.data.cpu().numpy(), loss2.data.cpu().numpy(), loss3.data.cpu().numpy()))


        loss.backward()

        #torch.nn.utils.clip_grad_norm(model.parameters(), 0.5)
        utils.ensure_shared_grads(model, shared_model, True)
        optimizer.step()

if __name__ == '__main__':

    parser = argparse.ArgumentParser() 
    parser.add_argument('--lr', help="learning rate", type=float, default = 1e-4)
    parser.add_argument('--num_epochs', help='number of epochs', type = int, default = 10000)
    parser.add_argument('--ncpus', help = 'number of cpus', type = int, default = 1) 
    parser.add_argument('--item_per_cycle', help = 'iteration per cycle', type = int, default = 128) 
    parser.add_argument('--dim_of_node_vector', help = 'dimension of node_vector', type = int, default = 128) 
    parser.add_argument('--dim_of_edge_vector', help = 'dimension of edge vector', type = int, default = 128) 
    parser.add_argument('--dim_of_FC', help = 'dimension of FC', type = int, default = 128) 
    parser.add_argument('--beta1', help = 'beta1: lambda paramter for VAE training', type = float, default = 5e-3) 
    parser.add_argument('--smiles_path', help='SMILES-data path')
    parser.add_argument('--data_paths', help='property-data paths', nargs='+', default=[], metavar='PATH')
    parser.add_argument('--save_dir', help = 'save directory', type = str) 
    parser.add_argument('--save_every', help = 'choose how often model will be saved', type = int, default = 200) 
    parser.add_argument('--shuffle_order', help = 'shuffle order or adding node and edge', action='store_true') 
    parser.add_argument('--active_ratio', help='active ratio in sampling (default: no matter)', type=float)
    parser.add_argument('--save_fpath', help='path of a saved model to restart') 
    args = parser.parse_args()

    # Process file/directory paths.
    depath = lambda path: os.path.realpath(os.path.expanduser(path))
    smiles_path = depath(args.smiles_path)
    data_paths = [depath(path) for path in args.data_paths]
    save_fpath = depath(args.save_fpath) if args.save_fpath else None
    save_dir = depath(args.save_dir)
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    # Load data.
    # See `utils.load_data` for the variable structures.
    id_to_smiles, id_to_whole_conditions, id_to_scaffold_conditions = utils.load_data(smiles_path, *data_paths)

    # Get the number of conditions as a hyperparameter.
    # `args` delivers the hyperparameters to `ggm.ggm`.
    args.N_conditions = len(next(iter(id_to_whole_conditions.values()))) + len(next(iter(id_to_scaffold_conditions.values())))

    #lines for multiprocessing
    mp.set_start_method('spawn')
    torch.manual_seed(1)

    #model 
    shared_model = ggm(args)
    shared_model.share_memory()  # torch.nn.Module.share_memory

    #shared optimizer
    shared_optimizer = SharedAdam(shared_model.parameters(), lr=args.lr, amsgrad=True)
    shared_optimizer.share_memory()
    print ("Model #Params: %dK" % (sum([x.nelement() for x in shared_model.parameters()]) / 1000,))
    
    #initialize parameters of the model 
    if save_fpath:
        initial_epoch, initial_cycle = [int(value) for value in re.findall('\d+', os.path.basename(save_fpath))]
        shared_model = utils.initialize_model(shared_model, save_fpath)
    else:
        initial_epoch = initial_cycle = 0
        shared_model = utils.initialize_model(shared_model, False)

    num_cycles = int(len(id_to_smiles)/args.ncpus/args.item_per_cycle)
    print(f"""\
ncpus             : {args.ncpus}
OMP_NUM_THREADS   : {os.environ.get('OMP_NUM_THREADS')}
Number of data    : {len(id_to_smiles)}
Number of epochs  : {args.num_epochs}
Number of cycles  : {num_cycles} per epoch
Minibatch size    : {args.item_per_cycle} per CPU per cycle
Learning rate     : {args.lr}
dim_of_node_vector: {args.dim_of_node_vector}
dim_of_edge_vector: {args.dim_of_edge_vector}
dim_of_FC         : {args.dim_of_FC}
beta1             : {args.beta1}
SMILES data path  : {smiles_path}
Data directories  : {data_paths}
Save directory    : {save_dir}
Save model every  : {args.save_every} cycles per epoch (Total {args.num_epochs*(num_cycles//args.save_every+1)} models)
shuffle_order     : {args.shuffle_order}
Restart from      : {save_fpath}
""")
    
    print("epoch  cyc  totcyc  loss  loss1  loss2  loss3  time")
    for epoch in range(args.num_epochs):
        # Jump to the previous epoch for restart.
        if epoch < initial_epoch:
            continue
        for cycle in range(num_cycles):
            # Jump to the previous cycle for restart.
            if epoch == initial_epoch:
                if cycle < initial_cycle:
                    continue
            retval_list = mp.Manager().list()  # Is this needed?
            # List of multiprocessing.managers.ListProxy to collect losses
            retval_list = [mp.Manager().list() for i in range(args.ncpus)]
            st = time.time()
            processes = []
            for pid in range(args.ncpus):
                # Sample keys without considering activeness.
                if args.active_ratio is None:
                    keys = random.sample(id_to_smiles.keys(), args.item_per_cycle)
                # Sample active and inactive keys by the required ratio.
                else:
                    keys = utils.sample_data(id_to_whole_conditions, args.item_per_cycle, args.active_ratio)

                # Property (descriptor) values work as conditions.
                # We need both of whole and scaffold values.
                # whole_conditions := [
                #     [ value1, value2, ... ],  # condition values of whole 1
                #     [ value1, value2, ... ],  # condition values of whole 2
                #     ... ]
                whole_conditions = [id_to_whole_conditions[key] for key in keys]
                scaffold_conditions = [id_to_scaffold_conditions[key] for key in keys]

                # SMILESs of whole molecules and scaffolds.
                wholes = [id_to_smiles[key][0] for key in keys]
                scaffolds = [id_to_smiles[key][1] for key in keys]

                proc = mp.Process(target=train, args=(shared_model, shared_optimizer, wholes, scaffolds, whole_conditions, scaffold_conditions, pid, retval_list, args))
                proc.start()
                processes.append(proc)
                time.sleep(0.1)
            for proc in processes:
                proc.join() 
            end = time.time()

            # retval_list shape -> (ncpus, item_per_cycle, 4),
            loss = np.mean(np.array([losses[0] for k in retval_list for losses in k]))
            loss1 = np.mean(np.array([losses[1] for k in retval_list for losses in k]))
            loss2 = np.mean(np.array([losses[2] for k in retval_list for losses in k]))
            loss3 = np.mean(np.array([losses[3] for k in retval_list for losses in k]))
            print ('%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' %(epoch, cycle, epoch*num_cycles+cycle, loss, loss1, loss2, loss3, end-st))
            if cycle%args.save_every == 0:
                name = save_dir+'/save_'+str(epoch)+'_' + str(cycle)+'.pt'
                torch.save(shared_model.state_dict(), name)
