import torch
import numpy as np
import abc
from tqdm import trange
import torch.nn.functional as F
from losses import get_score_fn
from utils.graph_utils import mask_adjs, mask_x, gen_noise
from sde import VPSDE, subVPSDE
import pdb

controller_loss = [[],[]]


class Predictor(abc.ABC):
  """The abstract class for a predictor algorithm."""
  def __init__(self, sde, score_fn, probability_flow=False):
    super().__init__()
    self.sde = sde
    # Compute the reverse SDE/ODE
    self.rsde = sde.reverse(score_fn, probability_flow)
    self.score_fn = score_fn

  @abc.abstractmethod
  def update_fn(self, x, t, flags):
    pass


class Corrector(abc.ABC):
  """The abstract class for a corrector algorithm."""
  def __init__(self, sde, score_fn, snr, scale_eps, n_steps):
    super().__init__()
    self.sde = sde
    self.score_fn = score_fn
    self.snr = snr
    self.scale_eps = scale_eps
    self.n_steps = n_steps

  @abc.abstractmethod
  def update_fn(self, x, t, flags):
    pass


class EulerMaruyamaPredictor(Predictor):
  def __init__(self, obj, sde, score_fn, probability_flow=False):
    super().__init__(sde, score_fn, probability_flow)
    self.obj = obj

  def update_fn(self, x, adj, flags, t):
    dt = -1. / self.rsde.N

    if self.obj=='x':
      z = gen_noise(x, flags, sym=False)
      drift, diffusion = self.rsde.sde(x, adj, flags, t, is_adj=False)
      x_mean = x + drift * dt
      x = x_mean + diffusion[:, None, None] * np.sqrt(-dt) * z
      return x, x_mean

    elif self.obj=='adj':
      z = gen_noise(adj, flags)
      drift, diffusion = self.rsde.sde(x, adj, flags, t, is_adj=True)
      adj_mean = adj + drift * dt
      adj = adj_mean + diffusion[:, None, None] * np.sqrt(-dt) * z

      return adj, adj_mean

    else:
      raise NotImplementedError(f"obj {self.obj} not yet supported.")


class ReverseDiffusionPredictor(Predictor):
  def __init__(self, obj, sde, score_fn, probability_flow=False):
    super().__init__(sde, score_fn, probability_flow)
    self.obj = obj

  def update_fn(self, x, adj, flags, t):

    if self.obj == 'x':
      f, G = self.rsde.discretize(x, adj, flags, t, is_adj=False)
      z = gen_noise(x, flags, sym=False)
      x_mean = x - f
      x = x_mean + G[:, None, None] * z
      return x, x_mean

    elif self.obj == 'adj':
      f, G = self.rsde.discretize(x, adj, flags, t, is_adj=True)
      z = gen_noise(adj, flags)
      adj_mean = adj - f
      adj = adj_mean + G[:, None, None] * z
      return adj, adj_mean
    
    else:
      raise NotImplementedError(f"obj {self.obj} not yet supported.")


class NoneCorrector(Corrector):
  """An empty corrector that does nothing."""

  def __init__(self, obj, sde, score_fn, snr, scale_eps, n_steps):
    self.obj = obj
    pass

  def update_fn(self, x, adj, flags, t):
    if self.obj == 'x':
      return x, x
    elif self.obj == 'adj':
      return adj, adj
    else:
      raise NotImplementedError(f"obj {self.obj} not yet supported.")


class LangevinCorrector(Corrector):
  def __init__(self, obj, sde, score_fn, snr, scale_eps, n_steps):
    super().__init__(sde, score_fn, snr, scale_eps, n_steps)
    self.obj = obj

  def update_fn(self, x, adj, flags, t):
    sde = self.sde
    score_fn = self.score_fn
    n_steps = self.n_steps
    target_snr = self.snr
    seps = self.scale_eps

    if isinstance(sde, VPSDE) or isinstance(sde, subVPSDE):
      timestep = (t * (sde.N - 1) / sde.T).long()
      alpha = sde.alphas.to(t.device)[timestep]
    else:
      alpha = torch.ones_like(t)

    if self.obj == 'x':
      for i in range(n_steps):
        
        grad = score_fn(x, adj, flags, t)
        noise = gen_noise(x, flags, sym=False)
        grad_norm = torch.norm(grad.reshape(grad.shape[0], -1), dim=-1).mean()
        noise_norm = torch.norm(noise.reshape(noise.shape[0], -1), dim=-1).mean()
        step_size = (target_snr * noise_norm / grad_norm) ** 2 * 2 * alpha
        x_mean = x + step_size[:, None, None] * grad
        x = x_mean + torch.sqrt(step_size * 2)[:, None, None] * noise * seps
      return x, x_mean

    elif self.obj == 'adj':
      for i in range(n_steps):
        grad = score_fn(x, adj, flags, t)
        noise = gen_noise(adj, flags)
        grad_norm = torch.norm(grad.reshape(grad.shape[0], -1), dim=-1).mean()
        noise_norm = torch.norm(noise.reshape(noise.shape[0], -1), dim=-1).mean()
        step_size = (target_snr * noise_norm / grad_norm) ** 2 * 2 * alpha
        adj_mean = adj + step_size[:, None, None] * grad
        adj = adj_mean + torch.sqrt(step_size * 2)[:, None, None] * noise * seps
      return adj, adj_mean

    else:
      raise NotImplementedError(f"obj {self.obj} not yet supported")


# -------- PC sampler --------
def get_pc_sampler(sde_x, sde_adj, shape_x, shape_adj, predictor='Euler', corrector='None', 
                   snr=0.1, scale_eps=1.0, n_steps=1, 
                   probability_flow=False, continuous=False,
                   denoise=True, eps=1e-3, device='cuda'):

  def pc_sampler(model_x, model_adj, init_flags):

    score_fn_x = get_score_fn(sde_x, model_x, train=False, continuous=continuous)
    score_fn_adj = get_score_fn(sde_adj, model_adj, train=False, continuous=continuous)

    predictor_fn = ReverseDiffusionPredictor if predictor=='Reverse' else EulerMaruyamaPredictor 
    corrector_fn = LangevinCorrector if corrector=='Langevin' else NoneCorrector

    predictor_obj_x = predictor_fn('x', sde_x, score_fn_x, probability_flow)
    corrector_obj_x = corrector_fn('x', sde_x, score_fn_x, snr, scale_eps, n_steps)

    predictor_obj_adj = predictor_fn('adj', sde_adj, score_fn_adj, probability_flow)
    corrector_obj_adj = corrector_fn('adj', sde_adj, score_fn_adj, snr, scale_eps, n_steps)

    with torch.no_grad():
      # -------- Initial sample --------
      x = sde_x.prior_sampling(shape_x).to(device) 
      adj = sde_adj.prior_sampling_sym(shape_adj).to(device) 
      flags = init_flags
      x = mask_x(x, flags)
      adj = mask_adjs(adj, flags)
      diff_steps = sde_adj.N
      timesteps = torch.linspace(sde_adj.T, eps, diff_steps, device=device)

      # -------- Reverse diffusion process --------
      for i in trange(0, (diff_steps), desc = '[Sampling]', position = 1, leave=False):
        t = timesteps[i]
        vec_t = torch.ones(shape_adj[0], device=t.device) * t

        _x = x
        x, x_mean = corrector_obj_x.update_fn(x, adj, flags, vec_t)
        adj, adj_mean = corrector_obj_adj.update_fn(_x, adj, flags, vec_t)

        _x = x
        x, x_mean = predictor_obj_x.update_fn(x, adj, flags, vec_t)
        adj, adj_mean = predictor_obj_adj.update_fn(_x, adj, flags, vec_t)
      print(' ')
      return (x_mean if denoise else x), (adj_mean if denoise else adj), diff_steps * (n_steps + 1)
  return pc_sampler


# -------- S4 solver --------
def S4_solver(sde_x, sde_adj, shape_x, shape_adj, predictor='None', corrector='None', 
                        snr=0.1, scale_eps=1.0, n_steps=1, 
                        probability_flow=False, continuous=False,
                        denoise=True, eps=1e-3, device='cuda'):

  def s4_solver(model_x, model_adj, init_flags):

    score_fn_x = get_score_fn(sde_x, model_x, train=False, continuous=continuous)
    score_fn_adj = get_score_fn(sde_adj, model_adj, train=False, continuous=continuous)

    with torch.no_grad():
      # -------- Initial sample --------
      x = sde_x.prior_sampling(shape_x).to(device) 
      adj = sde_adj.prior_sampling_sym(shape_adj).to(device) 
      flags = init_flags
      x = mask_x(x, flags)
      adj = mask_adjs(adj, flags)
      diff_steps = sde_adj.N
      timesteps = torch.linspace(sde_adj.T, eps, diff_steps, device=device)
      dt = -1. / diff_steps

      # -------- Rverse diffusion process --------
      for i in trange(0, (diff_steps), desc = '[Sampling]', position = 1, leave=False):
        t = timesteps[i]
        vec_t = torch.ones(shape_adj[0], device=t.device) * t
        vec_dt = torch.ones(shape_adj[0], device=t.device) * (dt/2) 

        # -------- Score computation --------
        score_x = score_fn_x(x, adj, flags, vec_t)
        score_adj = score_fn_adj(x, adj, flags, vec_t)

        Sdrift_x = -sde_x.sde(x, vec_t)[1][:, None, None] ** 2 * score_x
        Sdrift_adj  = -sde_adj.sde(adj, vec_t)[1][:, None, None] ** 2 * score_adj

        # -------- Correction step --------
        timestep = (vec_t * (sde_x.N - 1) / sde_x.T).long()

        noise = gen_noise(x, flags, sym=False)
        grad_norm = torch.norm(score_x.reshape(score_x.shape[0], -1), dim=-1).mean()
        noise_norm = torch.norm(noise.reshape(noise.shape[0], -1), dim=-1).mean()
        if isinstance(sde_x, VPSDE):
          alpha = sde_x.alphas.to(vec_t.device)[timestep]
        else:
          alpha = torch.ones_like(vec_t)
      
        step_size = (snr * noise_norm / grad_norm) ** 2 * 2 * alpha
        x_mean = x + step_size[:, None, None] * score_x
        x = x_mean + torch.sqrt(step_size * 2)[:, None, None] * noise * scale_eps

        noise = gen_noise(adj, flags)
        grad_norm = torch.norm(score_adj.reshape(score_adj.shape[0], -1), dim=-1).mean()
        noise_norm = torch.norm(noise.reshape(noise.shape[0], -1), dim=-1).mean()
        if isinstance(sde_adj, VPSDE):
          alpha = sde_adj.alphas.to(vec_t.device)[timestep] # VP
        else:
          alpha = torch.ones_like(vec_t) # VE
        step_size = (snr * noise_norm / grad_norm) ** 2 * 2 * alpha
        adj_mean = adj + step_size[:, None, None] * score_adj
        adj = adj_mean + torch.sqrt(step_size * 2)[:, None, None] * noise * scale_eps

        # -------- Prediction step --------
        x_mean = x
        adj_mean = adj
        mu_x, sigma_x = sde_x.transition(x, vec_t, vec_dt)
        mu_adj, sigma_adj = sde_adj.transition(adj, vec_t, vec_dt) 
        x = mu_x + sigma_x[:, None, None] * gen_noise(x, flags, sym=False)
        adj = mu_adj + sigma_adj[:, None, None] * gen_noise(adj, flags)
        
        x = x + Sdrift_x * dt
        adj = adj + Sdrift_adj * dt

        mu_x, sigma_x = sde_x.transition(x, vec_t + vec_dt, vec_dt) 
        mu_adj, sigma_adj = sde_adj.transition(adj, vec_t + vec_dt, vec_dt) 
        x = mu_x + sigma_x[:, None, None] * gen_noise(x, flags, sym=False)
        adj = mu_adj + sigma_adj[:, None, None] * gen_noise(adj, flags)

        x_mean = mu_x
        adj_mean = mu_adj
      print(' ')
      return (x_mean if denoise else x), (adj_mean if denoise else adj), 0
  return s4_solver



# # 定义一个损失函数，用于获取图片的特征，然后与提示文字的特征进行对比
def cldr_loss_fn(drug_cell_features, text_features, number_label_embed, number_fusion_pred_embed, ic50):
  
    # input_normed = torch.nn.functional.normalize(drug_cell_features.
    #    unsqueeze(1), dim=2)
    
    
    # embed_normed = torch.nn.functional.normalize(text_features.repeat(drug_cell_features.shape[0],1).
    #    unsqueeze(0), dim=2)

    # dists_text = (
    #     input_normed.sub(embed_normed).norm(dim=2).div(2).arcsin().pow(2).mul(2)
    # ).mean()
    
    
    input_normed_num = torch.nn.functional.normalize(number_fusion_pred_embed.
       unsqueeze(1), dim=2)
    
    # pdb.set_trace()
    number_label_embed = number_label_embed.unsqueeze(0).repeat(number_fusion_pred_embed.shape[0],1,1)
    pred = torch.bmm(number_label_embed, input_normed_num.permute(0,2,1))
    
    pred = F.softmax(pred.squeeze(), dim=1)
    
    # 
    dists_num =  F.cross_entropy(pred.cpu(), torch.full((pred.shape[0],pred.shape[1]),float(ic50/0.1)))
    # dists = dists_text * 0.2 + dists_num * 0.8
    # 使用Squared Great Circle Distance计算距离
    # dists_num = torch.cosine_similarity(input_normed_num, embed_normed_num, dim=2).mean()
   
    global controller_loss
    # controller_loss[0].append(dists_text.item())
    controller_loss[1].append(dists_num.item())
    if len(controller_loss[1])==4:
        print(str(np.mean(controller_loss[0]))[:6],'\t', str(np.mean(controller_loss[1]))[:6],'\t')
        controller_loss = [[],[]]
     
    # return torch.cosine_similarity(input_normed_num, embed_normed_num, dim=2).mean()
    return dists_num

# # 定义一个损失函数，用于获取图片的特征，然后与提示文字的特征进行对比
def mol_sim_loss_fn(drug_aim, drug_generated):

    # pdb.set_trace()

    # dists = dists_text * 0.2 + dists_num * 0.8
    # 使用Squared Great Circle Distance计算距离
    dists_num = torch.cosine_similarity(drug_aim.repeat(drug_generated.shape[0],1), drug_generated, dim=1).mean()
   
    global controller_loss
    controller_loss[0].append(dists_num.item())
    controller_loss[1].append(dists_num.item())
    if len(controller_loss[1])==4:
        print(str(np.mean(controller_loss[0]))[:6],'\t', str(np.mean(controller_loss[1]))[:6],'\t')
        controller_loss = [[],[]]
     
    # return torch.cosine_similarity(input_normed_num, embed_normed_num, dim=2).mean()
    return dists_num
  
  
  
# # 定义一个损失函数，用于获取图片的特征，然后与提示文字的特征进行对比
# def cldr_loss_fn(drug_cell_features, text_features):
#     DC_cross_entropy_loss = torch.nn.CrossEntropyLoss()
#     T_cross_entropy_loss = torch.nn.CrossEntropyLoss()
#     # input_normed = torch.nn.functional.normalize(drug_cell_features.
#     #    unsqueeze(1), dim=2)
    
#     # embed_normed = torch.nn.functional.normalize(text_features.
#     #    unsqueeze(0), dim=2)
    
#     # pdb.set_trace()
    
#     dists = (drug_cell_features@text_features.repeat(drug_cell_features.shape[0],1).T)
#     dists_T = dists.T
    
#     labels = torch.arange(drug_cell_features.shape[0]).long().to(dists.device) 
#     loss_dc = DC_cross_entropy_loss(dists, labels)
#     loss_t = T_cross_entropy_loss(dists_T, labels)
    
#     loss_CLIP = (loss_dc + loss_t)/2
    
#     global controller_loss
#     controller_loss.append(loss_CLIP.item())
#     if len(controller_loss)==4:
#         print(str(np.mean(controller_loss))[:6],'\t')
#         controller_loss = []
     
#     return loss_CLIP





# def get_controller_grad_fn(flag, controller, config):
#   """Create the gradient function for the controller in use of -conditional sampling. """
#   def grad_fn(x, adj, ve_noise_scale, labels, config):
#     def prob_fn(x, adj, ve_noise_scale, labels, config):
#       # from torch_geometric import data as DATA
#       # classifier_grad_fn(DATA.Data(x=torch.tensor([0,0,0])),[0],[0])
#       controller.to(x.device)

#       x = x.requires_grad_()
#       adj = adj.requires_grad_()
      
      
#       # pdb.set_trace()
#       aim_mol = torch.load("/home/lk/project/DALLE24Drug/MolGen/GDSS/data/409.pth")
      
#       from utils.data_loader_mol import MolDataset, get_transform_fn
#       import networkx as nx
      
#       aim_dataset = MolDataset([aim_mol], get_transform_fn('GDSCv2'))
#       # aim_mols_nx = [adj for x, adj, label in aim_dataset]
#       # pdb.set_trace()
#       # with torch.no_grad(), torch.cuda.amp.autocast():
#       drug_aim_embed = controller.get_drug_embed(torch.tensor([x for x, adj, label in aim_dataset][0]).unsqueeze(0).to(x.device), torch.tensor([adj for x, adj, label in aim_dataset][0]).unsqueeze(0).to(x.device), ve_noise_scale, False)
#       drug_pred_embed = controller.get_drug_embed(x, adj, ve_noise_scale, True)
#       # pdb.set_trace()
      

#       # prob = jax.nn.log_softmax(logits, axis=-1)[jnp.arange(labels.shape[0]), labels].sum()
#       return mol_sim_loss_fn(drug_aim_embed, drug_pred_embed)


#     if flag == "x":
#       # pdb.set_trace()
#       return torch.autograd.grad(outputs = prob_fn(x, adj, ve_noise_scale, labels, config), inputs = x)[0]
#     elif flag == "adj":
#       return torch.autograd.grad(outputs = prob_fn(x, adj, ve_noise_scale, labels, config), inputs = adj)[0]

#   return grad_fn

# def get_controller_grad_fn(flag, controller, config):
#   """Create the gradient function for the controller in use of -conditional sampling. """
#   def grad_fn(x, adj, ve_noise_scale, labels, config):
#     def prob_fn(x, adj, ve_noise_scale, labels, config):
#       # from torch_geometric import data as DATA
#       # classifier_grad_fn(DATA.Data(x=torch.tensor([0,0,0])),[0],[0])
#       controller.to(x.device)

#       x = x.requires_grad_()
#       adj = adj.requires_grad_()
      
#       # with torch.no_grad(), torch.cuda.amp.autocast():
#       text_label_embed, number_label_embed = controller.get_text_embed(x, adj, labels, config)
      
      
#       fusion_pred_embed, number_fusion_pred_embed = controller.get_fusion_embed(x, adj, ve_noise_scale, labels)

#       # prob = jax.nn.log_softmax(logits, axis=-1)[jnp.arange(labels.shape[0]), labels].sum()
#       return cldr_loss_fn(fusion_pred_embed, text_label_embed, number_label_embed, number_fusion_pred_embed, labels['ic50'])


#     if flag == "x":
#       # pdb.set_trace()
#       return torch.autograd.grad(outputs = prob_fn(x, adj, ve_noise_scale, labels, config), inputs = x)[0]
#     elif flag == "adj":
#       return torch.autograd.grad(outputs = prob_fn(x, adj, ve_noise_scale, labels, config), inputs = adj)[0]

#   return grad_fn

def get_controller_label_embed(controller, config_controller, batchsize, device):

  controller.eval()
  
  if config_controller['base_model'] == 'TransFGDRP':
    data= [config_controller['label']['frag']]
  elif config_controller['base_model'] == 'TransEDRP':
    data= [[torch.tensor([config_controller['label']['cell'],config_controller['label']['ic50']])]]
  elif config_controller['base_model'] == 'TransDTA':
    data= [[torch.tensor([config_controller['label']['target'],config_controller['label']['kd']])]]
  elif config_controller['base_model'] == 'TransDDI':
    data= [[[config_controller['label']['drug'],config_controller['label']['label']]]]

  # pdb.set_trace()
     
  condition_embeddings = controller.text_full_forward(data, device = device)
  no_condition_embeddings = controller.forward_null_text(data, device = device)
  # import pdb;pdb.set_trace()
  results = [ no_condition_embeddings.squeeze(2).repeat([batchsize,1]),  
                  condition_embeddings.squeeze(2).repeat([batchsize,1])
                ] 
  return results

# -------- PC condiftion sampler --------
def get_pc_conditional_sampler(controller, sde_x, sde_adj, shape_x, shape_adj, predictor='Euler', corrector='None', 
                   snr=0.1, scale_eps=1.0, n_steps=1, 
                   probability_flow=False, continuous=False,
                   denoise=True, eps=1e-3, config=None, device='cuda'):
  
  # classifier_grad_fn_x = get_controller_grad_fn('x', controller, config)
  # classifier_grad_fn_adj = get_controller_grad_fn('adj', controller, config)
  
  config_controller = config.controller
  labels = config.controller.label

  
  def pc_sampler(model_x, model_adj, init_flags, w):

    # score_fn_x = get_score_fn(sde_x, model_x, train=False, continuous=continuous)
    # score_fn_adj = get_score_fn(sde_adj, model_adj, train=False, continuous=continuous)

    # predictor_fn = ReverseDiffusionPredictor if predictor=='Reverse' else EulerMaruyamaPredictor 
    # corrector_fn = LangevinCorrector if corrector=='Langevin' else NoneCorrector

    # predictor_obj_x = predictor_fn('x', sde_x, score_fn_x, probability_flow)
    # corrector_obj_x = corrector_fn('x', sde_x, score_fn_x, snr, scale_eps, n_steps)

    # predictor_obj_adj = predictor_fn('adj', sde_adj, score_fn_adj, probability_flow)
    # corrector_obj_adj = corrector_fn('adj', sde_adj, score_fn_adj, snr, scale_eps, n_steps)


      
    
    score_fn_x = get_score_fn(sde_x, model_x, train=False, continuous=continuous)
    score_fn_adj = get_score_fn(sde_adj, model_adj, train=False, continuous=continuous)
    

    condition = get_controller_label_embed(controller, config_controller, shape_x[0], device)
    
    def conditional_corrector_update_fn_x(score_fn_x, x, adj, flags, vec_t):
      corrector_fn = LangevinCorrector if corrector=='Langevin' else NoneCorrector
      
      def total_grad_fn(x, adj, flags, t):
        # ve_noise_scale = sde_x.marginal_prob(x, t)[1]
       
        return score_fn_x(x, adj, condition[0], flags, t)*(1-w)+score_fn_x(x, adj, condition[1], flags, t)*(w)
      #

      corrector_obj_x = corrector_fn('x', sde_x, total_grad_fn, snr, scale_eps, n_steps)
      return corrector_obj_x.update_fn(x, adj, flags, vec_t)
      
    def conditional_corrector_update_fn_adj(score_fn_adj, x, adj, flags, vec_t):
      corrector_fn = LangevinCorrector if corrector=='Langevin' else NoneCorrector
      
      def total_grad_fn(x, adj, flags, t):
        # ve_noise_scale = sde_adj.marginal_prob(x, t)[1]
        return score_fn_adj(x, adj, condition[0], flags, t)*(1-w)+score_fn_adj(x, adj, condition[1], flags, t)*(w)
      #
      
      corrector_obj_adj = corrector_fn('adj', sde_adj, total_grad_fn, snr, scale_eps, n_steps)
      return corrector_obj_adj.update_fn(x, adj, flags, vec_t)
      
    def conditional_predictor_update_fn_x(score_fn_x, x, adj, flags, vec_t):
      predictor_fn = ReverseDiffusionPredictor if predictor=='Reverse' else EulerMaruyamaPredictor 
      
      def total_grad_fn(x, adj, flags, t):
        # ve_noise_scale = sde_x.marginal_prob(x, t)[1]
        # return score_fn_x(x, adj, flags, t) + 100*classifier_grad_fn_x(x, adj, ve_noise_scale, labels, config)
        return  score_fn_x(x, adj, condition[0], flags, t)*(1-w)+score_fn_x(x, adj, condition[1], flags, t)*(w)
       
      predictor_obj_x = predictor_fn('x', sde_x, total_grad_fn, probability_flow)
      return predictor_obj_x.update_fn(x, adj, flags, vec_t)
    
    def conditional_predictor_update_fn_adj(score_fn_adj, x, adj, flags, vec_t):
      predictor_fn = ReverseDiffusionPredictor if predictor=='Reverse' else EulerMaruyamaPredictor 
      
      def total_grad_fn(x, adj, flags, t):
        # ve_noise_scale = sde_adj.marginal_prob(x, t)[1]
        # return score_fn_adj(x, adj, flags, t) + 100*classifier_grad_fn_adj(x, adj, ve_noise_scale, labels, config)
        return score_fn_adj(x, adj, condition[0], flags, t)*(1-w)+score_fn_adj(x, adj, condition[1], flags, t)*(w)
      
      predictor_obj_adj = predictor_fn('adj', sde_adj, total_grad_fn, probability_flow)
      return predictor_obj_adj.update_fn(x, adj, flags, vec_t)
    
    
    with torch.no_grad():
        # -------- Initial sample --------
      x = sde_x.prior_sampling(shape_x).to(device) 
      adj = sde_adj.prior_sampling_sym(shape_adj).to(device) 
      flags = init_flags

      x = mask_x(x, flags)
      adj = mask_adjs(adj, flags)
      diff_steps = sde_adj.N
      # diff_steps = 100
      timesteps = torch.linspace(sde_adj.T, eps, diff_steps, device=device)

      # -------- Reverse diffusion process --------
      for i in trange(0, (diff_steps), desc = '[Sampling]', position = 1, leave=False):
        t = timesteps[i]
        vec_t = torch.ones(shape_adj[0], device=t.device) * t

        # _x = x
        # x, x_mean = corrector_obj_x.update_fn(x, adj, flags, vec_t)
        # adj, adj_mean = corrector_obj_adj.update_fn(_x, adj, flags, vec_t)

        # _x = x
        # x, x_mean = predictor_obj_x.update_fn(x, adj, flags, vec_t)
        # adj, adj_mean = predictor_obj_adj.update_fn(_x, adj, flags, vec_t)

        _x = x
        x, x_mean = conditional_corrector_update_fn_x(score_fn_x, x, adj, flags, vec_t)
        adj, adj_mean = conditional_corrector_update_fn_adj(score_fn_adj, _x, adj, flags, vec_t)
        
        _x = x
        x, x_mean = conditional_predictor_update_fn_x(score_fn_x, x, adj, flags, vec_t)
        adj, adj_mean = conditional_predictor_update_fn_adj(score_fn_adj, _x, adj, flags, vec_t)
        

        # import pdb;pdb.set_trace()
      #   if i%100==0:  
      #     print(x_mean.sum().item(),adj_mean.sum().item())
      # print(f'w is {w}.')
      return (x_mean if denoise else x), (adj_mean if denoise else adj), diff_steps * (n_steps + 1)
  return pc_sampler
