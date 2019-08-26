#This program runs KT-scalar several times changing parameters such as filename or number of points taken
import os

#modify data

default = "default"

variables = {"file_name":default, 
			"final_grid" : default,
			"start_grid" : default,
			"initial_time":default, 
			"final_time" : default,
			"data_steps" : default,
			"int_steps"  : default}


def set_variable(variable, value):
  variables[variable] = value

def set_filename(filename):
	set_variable("file_name", filename)
	
def set_start_grid(startgrid):
	set_variable("start_grid", "{}".format(startgrid))
def set_final_grid(finalgrid):
	set_variable("final_grid", "{}".format(finalgrid))

def set_initial_time(initime):
	set_variable("initial_time", "{}".format(initime))
def set_final_time(finaltime):
	set_variable("final_time", "{}".format(finaltime))
	
def set_data_steps(datasteps):
	set_variable("data_steps", "{}".format(datasteps))
def set_int_steps(intsteps):
	set_variable("int_steps", "{}".format(intsteps))	

	
def changedata(variables):
  data = open("data","r")
  lines = data.readlines()
  data.close()
  data = open("data", "w")
  default = "default"	
  for line in lines:
    for variable in variables.keys():
      if line.startswith(variable):
        if variables[variable] != default:
          line = variable+"="+variables[variable]+"\n"
    data.write(line)
  data.close()


def run():
  os.system("./ff_exec")
  
def make():
  os.system("make")
  
def define_periodic():
  first_macro = open("first_macro_1d.h", "r")
  lines = first_macro.readlines()
  first_macro.close()
  first_macro = open("first_macro_1d.h", "w")
  for line in lines:
    if line.startswith("#undef PERIODIC"):
      line = "//#undef PERIODIC\n"
    first_macro.write(line)
  first_macro.close()
  
def undefine_periodic():
  first_macro = open("first_macro_1d.h", "r")
  lines = first_macro.readlines()
  for line in lines:
    print(line)
  first_macro.close()
  first_macro = open("first_macro_1d.h", "w")
  for line in lines:
    if line.startswith("//#undef PERIODIC"):
      line = "#undef PERIODIC\n"
    first_macro.write(line)
  first_macro.close()
  
def make_interface():
  undefine_periodic()
  make()
  
def make_periodic():
  define_periodic()
  make()
