source config.bash
python3 brown_wnonbond.py -yamlFile "crowder_EK.yaml" -run
#python3 brown_wnonbond.py -yamlFile "crowder_EK2.yaml" -run

#Different number of crowder cases:
python3 brown_wnonbond.py -yamlFile crowder_EK_1_crowders.yaml -run # crowders=1
python3 brown_wnonbond.py -yamlFile crowder_EK_4_crowders.yaml -run # crowders=4
python3 brown_wnonbond.py -yamlFile crowder_EK_9_crowders.yaml -run # crowders=9
python3 brown_wnonbond.py -yamlFile crowder_EK_16_crowders.yaml -run # crowders=16
python3 brown_wnonbond.py -yamlFile crowder_EK_25_crowders.yaml -run # crowders=25
python3 brown_wnonbond.py -yamlFile crowder_EK_36_crowders.yaml -run # crowders=36
python3 brown_wnonbond.py -yamlFile crowder_EK_49_crowders.yaml -run # crowders=49


#Different crowder dimension cases:
python3 brown_wnonbond.py -yamlFile crowder_EK_25_crowderDim.yaml -run # crowders=16, crowderDim=25
python3 brown_wnonbond.py -yamlFile crowder_EK_50_crowderDim.yaml -run # crowders=16, crowderDim=50
python3 brown_wnonbond.py -yamlFile crowder_EK_75_crowderDim.yaml -run # crowders=16, crowderDim=75
python3 brown_wnonbond.py -yamlFile crowder_EK_125_crowderDim.yaml -run # crowders=16, crowderDim=125
python3 brown_wnonbond.py -yamlFile crowder_EK_150_crowderDim.yaml -run # crowders=16, crowderDim=150
python3 brown_wnonbond.py -yamlFile crowder_EK_175_crowderDim.yaml -run # crowders=16, crowderDim=175

#Different crowder radius cases:
python3 brown_wnonbond.py -yamlFile crowder_EK_25_crowderRad.yaml -run # crowders=16, crowderDim=175, crowderRad=25
python3 brown_wnonbond.py -yamlFile crowder_EK_50_crowderRad.yaml -run # crowders=16, crowderDim=175, crowderRad=50
python3 brown_wnonbond.py -yamlFile crowder_EK_75_crowderRad.yaml -run # crowders=16, crowderDim=175, crowderRad=75
