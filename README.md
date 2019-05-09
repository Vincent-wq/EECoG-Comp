# EECoG-Comp (organizing)
This is an open source matlab based tool box for simltaneous EEG and ECoG data analysis and simulation by Qing Wang et al.( Brain Topogr. paper).
The folders are organized as follows:
1. EECoG-FWM: the forward modeling codes (related external tools);
2. LF: the EEG and ECoG lead field output folder;
3. newSimBioModel: the input anatomical and electrode position files;
4. EECoG_PreProc: the preprocessing module;
5. Laplacian: the Laplacian calculation module;
6. tools: common tools used by this project;
The files ara organized as follows:
1. EECoG_Comp_Sim_min_test_Human.m: human simulation test;
2. EECoG_Comp_Sim_min_test_Monkey.m: monkey simulation test;
3. res_collectHumanSimulationBlock.m: collect block source human simulation results;
4. res_collectHumanSimulationRand.m: collect random source human simulation results;
5. res_collectMonkeySimulationBlock.m: collect block source Monkey simulation results;
6. res_collectMonkeySimulationRand.m: collect random source Monkey simulation results;
7. Framework.tif: the framework of this pipeline.