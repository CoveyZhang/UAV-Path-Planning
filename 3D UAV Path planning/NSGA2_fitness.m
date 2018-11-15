function [ fitness_fit_out ] = NSGA2_fitness( dna_fit_in )
%NSGA2_FITNESS Summary of this function goes here
dnanum_fit=size(dna_fit_in,1);
fitness_fit_out=zeros(dnanum_fit,3);
fitness_fit_out(:,1)=NSGA2_fitness1(dna_fit_in);
fitness_fit_out(:,2)=NSGA2_fitness2(dna_fit_in);
fitness_fit_out(:,3)=NSGA2_fitness3(dna_fit_in);
%fitness_fit_out(:,3)=zeros(dnanum_fit,1);
%   Detailed explanation goes here
