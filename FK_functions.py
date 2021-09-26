import numpy as np
from numpy import sin, cos, pi

def r_z(q):
  rot = np.zeros((4,4),dtype=np.float64)
  rot[3,3] = 1
  rot[2,2] = 1
  rot[0,:] = [cos(q), -sin(q), 0, 0]
  rot[1,:] = [sin(q), cos(q), 0, 0]

  return rot

def r_y(q):
  rot = np.zeros((4,4),dtype=np.float64)
  rot[3,3] = 1
  rot[1,1] = 1
  rot[0,:] = [cos(q), 0, sin(q), 0]
  rot[2,:] = [-sin(q), 0, cos(q), 0]

  return rot

def r_x(q):
  rot = np.zeros((4,4),dtype=np.float64)
  rot[3,3] = 1
  rot[0,0] = 1
  rot[1,:] = [0, cos(q), -sin(q), 0]
  rot[2,:] = [0, sin(q), cos(q), 0]

  return rot

def tr(vector):
  mat = np.zeros((4,4),dtype=np.float64)
  mat[0:3,3] = vector
  np.fill_diagonal(mat,1)  

  return mat

def frame_trform(q_params): 

  z, y, x, t_z, t_y, t_x = q_params
  frame  = np.linalg.multi_dot([tr([0,0,t_z]), r_z(z), tr([0,t_y,0]), r_y(y), tr([t_x,0,0]), r_x(x) ])
  return frame

def serial_transform(Q):

  j = np.identity(4)

  frames = [j]

  for q_params in Q:
    j = np.dot(j, frame_trform(q_params))
    frames.append(j)
  
  return frames

def make_table(q_angles):
  q_1, q_2, q_3, q_4, q_5, q_6 = q_angles

  t_base = 0
  t_d0 = 0
  t_link_1 = 0.5
  t_link_2 = 0.5
  t_link_3 = 0.5
  t_link_4 = 0
  t_link_6 = 0.1
  t_tool = 0.1


  frame_1 = [q_1, 0, 0, t_base + t_d0 +t_link_1, 0, 0]
  frame_2 = [0, q_2, 0, 0, 0, t_link_2]
  frame_3 = [0, q_3, 0, 0, 0, t_link_3]
  frame_4 = [0, 0, q_4, t_link_4, 0, 0]
  frame_5 = [0, q_5, 0, 0, 0, 0]
  frame_6 = [0, 0, q_6, 0, 0, t_link_6 + t_tool]

  frame_0_to_6 = [frame_1, frame_2, frame_3, frame_4, frame_5, frame_6]
  return frame_0_to_6






