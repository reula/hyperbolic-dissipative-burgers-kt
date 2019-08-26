from pythonrunner_functions import *

set_filename("adv_cnvg_per")
set_start_grid(0)
set_initial_time(0.0)
set_final_time(2.0)
set_data_steps(20)
set_int_steps(200)
set_final_grid(400)

finalgridsperiodic=[200,400,800,1600]

make_periodic()

for finalgrid in finalgridsperiodic:
  set_filename("adv_cnvg_per")
  set_final_grid(finalgrid)
  changedata(variables)
  run()
