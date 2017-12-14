function make(arg1) 
% This make.m is used for Matlab under Windows

%	$Id: make.m 8223 2011-04-20 09:16:06Z schloegl $
%	Copyright (C) 2010,2011 by Alois Schloegl <alois.schloegl@gmail.com>
%       This function is part of the NaN-toolbox
%       http://pub.ist.ac.at/~schloegl/matlab/NaN/

% add -largeArrayDims on 64-bit machines

if (nargin>0 && strcmp(arg1,'clean')),
	if strcmp(computer,'PCWIN')
		dos('del *.obj');
		dos('del *.mex*');
	else
		unix('rm *.o'); 
		unix('rm *.mex*');
	end;
	return; 
end;

mex covm_mex.cpp
mex sumskipnan_mex.cpp
mex histo_mex.cpp
mex kth_element.cpp
mex str2array.cpp
mex xptopen.cpp
mex -c svm.cpp
mex -c svm_model_matlab.c
mex -c tron.cpp
mex -c linear.cpp
mex -c linear_model_matlab.c
if strcmp(computer,'PCWIN') && ~exist('OCTAVE_VERSION','builtin'),
	mex svmtrain_mex.cpp svm.obj svm_model_matlab.obj
	mex svmpredict_mex.cpp svm.obj svm_model_matlab.obj

        if ~exist('LAPACK/daxpy.f','file') || ~exist('LAPACK/ddot.f','file') || ~exist('LAPACK/dscal.f','file') || ~exist('LAPACK/dnrm2.f','file'),
                fprintf(1,'The lapack functions daxpy, ddot, dscal, and dnrm2 are required.\n');
                fprintf(1,'If some functions are missing, get them from here:\n');
                if ~exist('LAPACK','dir') mkdir('LAPACK'); end; 
	        fprintf(1,'Get http://www.netlib.org/blas/daxpy.f and save to %s',fullfile(pwd,'LAPACK')); 
	        fprintf(1,'Get http://www.netlib.org/blas/ddot.f and save to %s',fullfile(pwd,'LAPACK')); 
	        fprintf(1,'Get http://www.netlib.org/blas/dscal.f and save to %s',fullfile(pwd,'LAPACK')); 
	        fprintf(1,'Get http://www.netlib.org/blas/dnrm2.f and save to %s',fullfile(pwd,'LAPACK')); 
                fprintf(1,'Press any key to continue ... '\n);
	        pause;
	end;         
        mex -c LAPACK/daxpy.f
        mex -c LAPACK/ddot.f
        mex -c LAPACK/dscal.f
        mex -c LAPACK/dnrm2.f
	dos('copy train.c train.cpp');
	mex('train.cpp','tron.obj','linear.obj','linear_model_matlab.obj','daxpy.obj','ddot.obj','dscal.obj','dnrm2.obj')
	dos('del *.obj');

else
	mex svmtrain_mex.cpp svm.o svm_model_matlab.o 
	mex svmpredict_mex.cpp svm.o svm_model_matlab.o
	unix('cp train.c train.cpp');
	mex train.cpp tron.o linear.o linear_model_matlab.o 
	unix('rm *.o');
end

