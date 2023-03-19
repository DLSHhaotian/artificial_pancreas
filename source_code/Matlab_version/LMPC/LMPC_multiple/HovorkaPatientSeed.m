function [seed,num] = HovorkaPatientSeed(numPatient)
% Syntax: [seed,num] = HovorkaPatientSeed(numPatient)
%         numPatient: number of patient
seed_vec = [18;27;39;42;72;420;612;745;775;866];
num_vec = [4;5;4;5;3;5;4;3;6;4];
seed = seed_vec(numPatient);
num = num_vec(numPatient);
end