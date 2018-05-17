%% Regression_test script to run slp_trust and sqp examples

%% Log file
clear; clc; close all
delete regression_test.txt
diary  regression_test.txt

%% Run scripts
run2barFox
runBeamGVslp
runHaftka4p2p1
runHaftka6p3p1slp
runRosenSuzuki
runSpringExampleTrust
runSvanbergSLP
runSvanbergSQP
runBeamGV
tsuite

diary off