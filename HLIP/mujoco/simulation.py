import mujoco as mj
import numpy as np
import os
from mujoco.glfw import glfw

# define the XMl file path
xml_file_path = "./models/achilles.xml"
absolute_xml_file_path = os.path.abspath(xml_file_path)

# create model from the XML file
model = mj.MjModel.from_xml_path(absolute_xml_file_path)

# print out some mujoco data
print("Model name: ", model.names)
print("Model nq: ", model.nq)
print("Model nv: ", model.nv)
print("Model nu: ", model.nu)
print("Model nbody: ", model.nbody)
print("Model njnt: ", model.njnt)
print("Model ngeom: ", model.ngeom)




