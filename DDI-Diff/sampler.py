import os
import time
import pickle
import math
import torch

from utils.logger import Logger, set_log, start_log, train_log, sample_log, check_log
from utils.loader import load_ckpt, load_data, load_seed, load_device, load_model_from_ckpt, load_model_params, \
                         load_ema_from_ckpt, load_sampling_fn, load_condition_sampling_fn, load_eval_settings
from utils.graph_utils import adjs_to_graphs, init_flags, quantize, quantize_mol
from utils.plot import save_graph_list, plot_graphs_list
from evaluation.stats import eval_graph_list
from utils.mol_utils import gen_mol, mols_to_smiles, load_smiles, canonicalize_smiles, mols_to_nx, filter_smiles_with_labels

import sys
sys.path.insert(0,'/home/hch/DDI-Diff/DDI/moses')
from moses.metrics.metrics import get_all_metrics
from utils.mol_utils import mols_to_nx, smiles_to_mols
import pdb
from tqdm import tqdm 


# -------- Sampler for generic graph generation tasks --------
class Sampler(object):
    def __init__(self, config):
        super(Sampler, self).__init__()

        self.config = config
        self.device = load_device()

    def sample(self):
        # -------- Load checkpoint --------
        self.ckpt_dict = load_ckpt(self.config, self.device)
        self.configt = self.ckpt_dict['config']

        load_seed(self.configt.seed)

        self.train_graph_list, self.test_graph_list = load_data(self.configt, get_graph_list=True)

        self.log_folder_name, self.log_dir, _ = set_log(self.configt, is_train=False)
        self.log_name = f"{self.config.ckpt}-sample"
        logger = Logger(str(os.path.join(self.log_dir, f'{self.log_name}.log')), mode='a')

        if not check_log(self.log_folder_name, self.log_name):
            logger.log(f'{self.log_name}')
            start_log(logger, self.configt)
            train_log(logger, self.configt)
        sample_log(logger, self.config)

        # -------- Load models --------
        self.model_x = load_model_from_ckpt(self.ckpt_dict['params_x'], self.ckpt_dict['x_state_dict'], self.device)
        self.model_adj = load_model_from_ckpt(self.ckpt_dict['params_adj'], self.ckpt_dict['adj_state_dict'], self.device)

        if self.config.sample.use_ema:
            self.ema_x = load_ema_from_ckpt(self.model_x, self.ckpt_dict['ema_x'], self.configt.train.ema)
            self.ema_adj = load_ema_from_ckpt(self.model_adj, self.ckpt_dict['ema_adj'], self.configt.train.ema)
            
            self.ema_x.copy_to(self.model_x.parameters())
            self.ema_adj.copy_to(self.model_adj.parameters())

        self.sampling_fn = load_sampling_fn(self.configt, self.config.sampler, self.config.sample, self.device)

        # -------- Generate samples --------
        logger.log(f'GEN SEED: {self.config.sample.seed}')
        load_seed(self.config.sample.seed)

        num_sampling_rounds = math.ceil(len(self.test_graph_list) / self.configt.data.batch_size)
        gen_graph_list = []
        for r in range(num_sampling_rounds):
            t_start = time.time()

            self.init_flags = init_flags(self.train_graph_list, self.configt).to(self.device[0])

            x, adj, _ = self.sampling_fn(self.model_x, self.model_adj, self.init_flags)

            logger.log(f"Round {r} : {time.time()-t_start:.2f}s")

            samples_int = quantize(adj)
            gen_graph_list.extend(adjs_to_graphs(samples_int, True))

        gen_graph_list = gen_graph_list[:len(self.test_graph_list)]

        # -------- Evaluation --------
        methods, kernels = load_eval_settings(self.config.data.data)
        result_dict = eval_graph_list(self.test_graph_list, gen_graph_list, methods=methods, kernels=kernels)
        logger.log(f'MMD_full {result_dict}', verbose=False)
        logger.log('='*100)

        # -------- Save samples --------
        save_dir = save_graph_list(self.log_folder_name, self.log_name, gen_graph_list)
        with open(save_dir, 'rb') as f:
            sample_graph_list = pickle.load(f)
        plot_graphs_list(graphs=sample_graph_list, title=f'{self.config.ckpt}', max_num=16, save_dir=self.log_folder_name)


# -------- Sampler for molecule generation tasks --------
class Sampler_mol(object):
    def __init__(self, config):
        self.config = config
        self.device = load_device()

    def sample(self):
        # -------- Load checkpoint --------
        self.ckpt_dict = load_ckpt(self.config, self.device)
        self.configt = self.ckpt_dict['config']

        load_seed(self.config.seed)

        self.log_folder_name, self.log_dir, _ = set_log(self.configt, is_train=False)
        self.log_name = f"{self.config.ckpt}-sample"
        logger = Logger(str(os.path.join(self.log_dir, f'{self.log_name}.log')), mode='a')

        if not check_log(self.log_folder_name, self.log_name):
            start_log(logger, self.configt)
            train_log(logger, self.configt)
        sample_log(logger, self.config)

        # -------- Load models --------
        self.model_x = load_model_from_ckpt(self.ckpt_dict['params_x'], self.ckpt_dict['x_state_dict'], self.device)
        self.model_adj = load_model_from_ckpt(self.ckpt_dict['params_adj'], self.ckpt_dict['adj_state_dict'], self.device)
        
        self.sampling_fn = load_sampling_fn(self.configt, self.config.sampler, self.config.sample, self.device)

        # -------- Generate samples --------
        logger.log(f'GEN SEED: {self.config.sample.seed}')
        load_seed(self.config.sample.seed)

        train_smiles, test_smiles = load_smiles(self.configt.data.data)
        train_smiles, test_smiles = canonicalize_smiles(train_smiles), canonicalize_smiles(test_smiles)

        self.train_graph_list, _ = load_data(self.configt, get_graph_list=True)     # for init_flags
        with open(f'{self.configt.data.dir}/{self.configt.data.data.lower()}_test_nx.pkl', 'rb') as f:
            self.test_graph_list = pickle.load(f)                                   # for NSPDK MMD
        # self.init_flags = init_flags(self.train_graph_list, self.configt, 1000).to(self.device[0])
        self.init_flags = init_flags(self.train_graph_list, self.configt, 1000).to(self.device[0])
        x, adj, _ = self.sampling_fn(self.model_x, self.model_adj, self.init_flags)
        
        samples_int = quantize_mol(adj)

        samples_int = samples_int - 1
        samples_int[samples_int == -1] = 3      # 0, 1, 2, 3 (no, S, D, T) -> 3, 0, 1, 2

        adj = torch.nn.functional.one_hot(torch.tensor(samples_int), num_classes=4).permute(0, 3, 1, 2)
        x = torch.where(x > 0.5, 1, 0)
        x = torch.concat([x, 1 - x.sum(dim=-1, keepdim=True)], dim=-1)      # 32, 9, 4 -> 32, 9, 5

        gen_mols, num_mols_wo_correction = gen_mol(x, adj, self.configt.data.data)
        num_mols = len(gen_mols)

        gen_smiles = mols_to_smiles(gen_mols)
        gen_smiles = [smi for smi in gen_smiles if len(smi)]
        
        # -------- Save generated molecules --------
        with open(os.path.join(self.log_dir, f'{self.log_name}.txt'), 'a') as f:
            for smiles in gen_smiles:
                f.write(f'{smiles}\n')

        # -------- Evaluation --------
        scores = get_all_metrics(gen=gen_smiles, k=len(gen_smiles), device=self.device[0], n_jobs=8, test=test_smiles, train=train_smiles)
        scores_nspdk = eval_graph_list(self.test_graph_list, mols_to_nx(gen_mols), methods=['nspdk'])['nspdk']

        logger.log(f'Number of molecules: {num_mols}')
        logger.log(f'validity w/o correction: {num_mols_wo_correction / num_mols}')
        for metric in ['valid', f'unique@{len(gen_smiles)}', 'FCD/Test', 'Novelty']:
            logger.log(f'{metric}: {scores[metric]}')
        logger.log(f'NSPDK MMD: {scores_nspdk}')
        logger.log('='*100)




# -------- Sampler for molecule generation tasks --------
class Sampler_mol_condition_clddi(object):
    def __init__(self, config, w=None):
        self.config = config
        self.device = load_device()
        self.params_x, self.params_adj = load_model_params(self.config)
        self.samples_num = 10
        self.w = 0.0 if w is None else w
        print("self.w is ", self.w)
        
    def sample(self):
        # -------- Load checkpoint --------
        self.ckpt_dict = load_ckpt(self.config, self.device)
        # self.ckpt_dict_condition = load_ckpt(self.config, self.device, market='')
        self.configt = self.ckpt_dict['config']
    
        load_seed(self.config.seed)

        self.log_folder_name, self.log_dir, _ = set_log(self.configt, is_train=False)
        self.log_name = f"{self.config.ckpt}-sample-{str(self.config.controller.label['drug'])}-{str(self.config.controller.label['label'])}"

        logger = Logger(str(os.path.join(self.log_dir, f'{self.log_name}.log')), mode='a')

        if not check_log(self.log_folder_name, self.log_name):
            start_log(logger, self.configt)
            train_log(logger, self.configt)
        sample_log(logger, self.config)


       
        
        # -------- Load models --------
        self.model_x = load_model_from_ckpt(self.ckpt_dict['params_x'], self.ckpt_dict['x_state_dict'], self.device, config_train=self.configt.train)
        self.model_adj = load_model_from_ckpt(self.ckpt_dict['params_adj'], self.ckpt_dict['adj_state_dict'], self.device, config_train=self.configt.train)
        
        # self.model_x_condition = load_model_from_ckpt(self.ckpt_dict_condition['params_x'], self.ckpt_dict_condition['x_state_dict'], self.device)
        # self.model_adj_condition = load_model_from_ckpt(self.ckpt_dict_condition['params_adj'], self.ckpt_dict_condition['adj_state_dict'], self.device)

        # 这里加循环 多生成几类条件对应的分子
        import pandas as pd
        drugB_and_label = pd.read_csv('/home/hch/DDI-Diff/DDI/data/filtered_clean_max_smiles.csv')
        matched_ddi = pd.read_csv('/home/hch/DDI-Diff/DDI/data/matched_ddi.csv')
        # self.log_dir = '/home/hch/DDI-Diff/DDI/results/'
        # self.log_dir = '/home/hch/DDI-Diff/DDI/results_drugbank_qm9/'
        # self.log_dir = '/home/hch/DDI-Diff/DDI/results_drugbank_qm9_35/'
        # self.log_dir = '/home/hch/DDI-Diff/DDI/results_drugbank_qm9_45/'
        self.log_dir = '/home/hch/DDI-Diff/DDI/results_w5.0/'
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        # 已经采了前20类（0--29）条件对应的分子，再往后采10类（30--38）
        # for i in tqdm(range(30, 39), desc="Generating molecules"):

        for i in tqdm([2, 26, 5, 23], desc="Generating molecules"):
            drug  = drugB_and_label.iloc[i]['Drug2']
            label = drugB_and_label.iloc[i]['Label']

            try:
                ddi_info = matched_ddi[matched_ddi['Interaction type'] == label]['DDI type ID'].iloc[0]
            except IndexError:
                print(f"No matched DDI type ID for label {label}")
                continue

            self.config['controller']['label']['drug'] = drug
            self.config['controller']['label']['label'] = label

            self.sampling_fn = load_condition_sampling_fn(self.configt, self.config, self.config.sampler, self.config.sample, self.device, self.params_x, self.params_adj, self.samples_num)

            # -------- Generate samples --------
            logger.log(f'GEN SEED: {self.config.sample.seed}')
            load_seed(self.config.sample.seed)

            train_smiles, _ = load_smiles(self.configt.data.data)
        
            test_topK_df_1 = filter_smiles_with_labels(self.config, topk=-1)
            # denote -1 as no number limitation.

            
            train_smiles = canonicalize_smiles(train_smiles)
        

            test_smiles_1 = canonicalize_smiles(test_topK_df_1)

        
            # self.train_graph_list, _ = load_data(self.configt, get_graph_list=True)     # for init_flags
            # with open(f'{self.configt.data.dir}/{self.configt.data.data.lower()}_test_nx.pkl', 'rb') as f:
            #     self.test_graph_list = pickle.load(f)                                   # for NSPDK MMD

            # self.init_flags = init_flags(self.train_graph_list, self.configt, self.samples_num).to(self.device[0])
            self.init_flags = torch.load("./temp/temp_data/init_flags_1000_new.pth")
            # torch.save(self.init_flags, "./temp/temp_data/init_flags_1000_new.pth")
            # import pdb;pdb.set_trace()
            # sys.exit()  # 程序在此结束
            # Deal with the self.test_graph_list as test_smiles(test_topK_df)

            self.test_topK_df_nx_graphs_1 = mols_to_nx(smiles_to_mols(test_smiles_1))
            # pdb.set_trace()
        
            x, adj, _ = self.sampling_fn(self.model_x, self.model_adj, self.init_flags, self.w)
            # x, adj, _ = self.sampling_fn(self.model_x, self.model_adj, None)
            # continue
            samples_int = quantize_mol(adj)

            samples_int = samples_int - 1
            samples_int[samples_int == -1] = 3      # 0, 1, 2, 3 (no, S, D, T) -> 3, 0, 1, 2
            # import pdb;pdb.set_trace()
            adj = torch.nn.functional.one_hot(torch.tensor(samples_int), num_classes=4).permute(0, 3, 1, 2)
            x = torch.where(x > 0.5, 1, 0)
            x = torch.concat([x, 1 - x.sum(dim=-1, keepdim=True)], dim=-1)      # 32, 9, 4 -> 32, 9, 5

            gen_mols, num_mols_wo_correction = gen_mol(x, adj, self.configt.data.data[0] if  type(self.configt.data.data) == list else self.configt.data.data)
            num_mols = len(gen_mols)

            gen_smiles = mols_to_smiles(gen_mols)
            gen_smiles = [smi for smi in gen_smiles if len(smi)]

            # -------- Save generated molecules --------
            # -------- 去重并创建 DataFrame --------
            unique_gen_smiles = list(set(gen_smiles))  # 去重
            drug_a_b_data = pd.DataFrame({
                'drug_A': unique_gen_smiles,
                'drug_B': [drug] * len(unique_gen_smiles)  # 复制 drug 列
            })

            ddi_info = matched_ddi[matched_ddi['Interaction type'] == label]['DDI type ID'].iloc[0]

            file_name = f'drug{i}?{label}?{ddi_info}.csv'
            drug_a_b_data.to_csv(os.path.join(self.log_dir, file_name), index=False)
            
         
            # -------- Evaluation --------
            try:
                scores_1 = get_all_metrics(gen=unique_gen_smiles, k=len(unique_gen_smiles), device=self.device[0], n_jobs=8, test=test_smiles_1, train=train_smiles)
                # pdb.set_trace()
                scores_nspdk_1 = eval_graph_list(self.test_topK_df_nx_graphs_1, mols_to_nx(smiles_to_mols(unique_gen_smiles)), methods=['nspdk'])['nspdk']

                with open(os.path.join(self.log_dir, f'drug{i}?{label}?{ddi_info}.txt'), 'w') as f:
                    f.write(f'drug_B: {drug}\n')
                    f.write(f'Interaction type: {label}\n')
                    f.write(f'DDI type ID: {ddi_info}\n')
                    f.write(f'Number of molecules: {num_mols}\n')
                    f.write(f'validity w/o correction: {num_mols_wo_correction / num_mols}\n')
                    f.write(f'==================: scores_1 ==================\n')
                    for metric in scores_1.keys():
                        f.write(f'{metric}: {scores_1[metric]}\n')
                    f.write(f'==================: scores_1 ==================\n')
                    f.write(f'NSPDK MMD: {scores_nspdk_1}\n')
                    f.write(f'number of test smiles is: {len(test_smiles_1)}\n')
                    f.write('='*100)
            except Exception as e:
                # try:
                #     os.remove(os.path.join(self.log_dir, file_name))
                # except FileNotFoundError:
                #     pass

                try:
                    os.remove(os.path.join(self.log_dir, f'drug{i}?{label}?{ddi_info}.txt'))
                except FileNotFoundError:
                    pass

                print(f"Error occurred: {e}")

        
