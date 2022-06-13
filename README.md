Monte Carlo Dose Calculation for Interventional Radiology (McDCIR)  
>> Dose Calculation Mode: Post-procedure Dose Calculation (PDC)
Contributor: Haegin Han, Sungho Moon, Gahee Son    

* Build (The execution file will be created in CMakeLists.txt path)
  1. mkdir build && cd build
  2. cmake ..
  3. make

* Configuration file
  1. PDC.dat:    Web-server connection information file (binary).
  2. PDC.server: Set server execution parameters such as server id, thread #, and nps/frame.
  3. PDC.column: Configuration of indices of DB table columnes in web-server. 
   'Radiologist_wxyz', 'DoseResults' should be configured like '(start column idx)-(end column idx)'.

* Phantom
  1. .dose: Do not revise the organ names and orders.

* Execution
./calculator [option 1] [option 2] [option 3] ...  
[options]
-m [G4 macro]         : Geant4 run macro file
-o [Output]           : Dose results file name
-p [Patient phantom]  : Patient phantom name (If patient phantom is not set up, patient will be not installed)
-d [Doctor phantom]   : Doctor phantom name (If doctor phantom is not set up, doctor will be not installed)  
-v                    : Qt Visualization (Don't use it with '-m' option)  
-sql                  : Run PDC system (This option will be default mode, Don't use it with '-m' option)
-id [Server id]       : Read server # parameter of 'PDC.server' and set them in PDC system. (Use it with '-sql' option)  
-x                    : Reset execution parameters in web-server (for testing)

example: ./calculator -m run.mac -o output -p RANDO_M -d MRCP_AM 
example: ./calculator -p RANDO_M -d MRCP_AM -v
example: ./calculator -p RANDO_M -d MRCP_AM -id 1 -sql
