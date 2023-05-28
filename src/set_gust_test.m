function gust_tests_data = set_gust_test(gust_test_case)

% Set base filepath
base_filepath = "../Gust Data/";

% Get gust test conditions according to case identifier
switch gust_test_case
    case 1
        filename = "1CG_H1PI_A0.mat";
    case 2
        filename = "1CG_H1PI_A10.mat";
    case 3
        filename = "1CG_H1PI_A15.mat";
    case 4
        filename = "1CG_H8PI_A0.mat";
    case 5
        filename = "1CG_H8PI_A10.mat";
    case 6
        filename = "1CG_H8PI_A15.mat";
    case 7
        filename = "SEG_A00.mat";
    case 8
        filename = "SEG_A10.mat";
    case 9
        filename = "SEG_A15.mat";
    case 10
        filename = "SG_A00.mat";
    case 11
        filename = "1CG_H_37_A09.mat";
    case 12
        filename = "1CG_H_45_A09.mat";
    case 13
        filename = "1CG_H_15_A11.mat";
    case 14
        filename = "SEG_M05_lambda05.mat";
    case 15
        filename = "SEG_M05_lambda06.mat";
    case 16
        filename = "SEG_M05_lambda07.mat";
    case 17
        filename = "SEG_M05_lambda075.mat";
    case 18
        filename = "SEG_M05_lambda1.mat";
    case 19
        filename = "SEG_M05_lambda15.mat";
    case 20
        filename = "SEG_M05_lambda2.mat";
    case 21
        filename = "SEG_inc_lambda0.mat"; 
    case 22
        filename = "SEG_inc_lambda025.mat"; 
    case 23
        filename = "SEG_inc_lambda05.mat";
    case 24
        filename = "SEG_inc_lambda075.mat"; 
    case 25
        filename = "SEG_inc_lambda1.mat"; 
    case 26
        filename = "SEG_inc_lambda15.mat";  
    case 27
        filename = "SEG_inc_lambda2.mat"; 
    case 28
        filename = "SEG_inc_lambdan025.mat";   
    case 29
        filename = "SEG_inc_lambdan05.mat";     
    case 30
        filename = "SEG_inc_lambdan075.mat"; 
    case 31
        filename = "SEG_inc_lambdan1.mat";
    case 32
        filename = "SEG_inc_lambdan15.mat";  
    case 33
        filename = "SEG_inc_lambdan2.mat";    
end

% Load data from current case and set on output cell
filepath = base_filepath + filename;
load(filepath,'gust_tests');
gust_tests_data = gust_tests{1};

end