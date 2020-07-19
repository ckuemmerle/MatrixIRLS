clear; close all; clc

% ========= random tests =======
fprintf('--- random matrix completion ---\n')

% --------- problem size -------
m = 1000;
n = 1000;
p = round(0.1*m*n);
r = 20;

fprintf('problem size: m = %d, n = %d, p = %d, r = %d\n\n', m, n, p, r)

% --------- NIHT_Matrix ---------
fprintf('NIHT starts ...\n')
matrix_completion('NIHT_Matrix', m, n, p, r)
fprintf('NIHT over ...\n\n')

% ---------- CGIHT_Matrix -------
fprintf('CGIHT starts ...\n')
matrix_completion('CGIHT_Matrix', m, n, p, r)
fprintf('CGIHT over ...\n\n')

% ----------- ASD --------------
fprintf('ASD starts ...\n')
matrix_completion('ASD', m, n, p, r)
fprintf('ASD over ...\n\n')

% -------- ScaledASD ------------
fprintf('ScaledASD starts ...\n')
matrix_completion('ScaledASD', m, n, p, r)
fprintf('ScaledASD over ...\n\n')

% ========= image inpainting ======
fprintf('--- random/cross inpainting ---\n')

alg_list = cell(4,1);
alg_list{1} = 'NIHT_Matrix';
alg_list{2} = 'CGIHT_Matrix';
alg_list{3} = 'ASD';
alg_list{4} = 'ScaledASD';

img_list = cell(3,1);
img_list{1} = 'boat.png';
img_list{2} = 'barb.png';
img_list{3} = 'lena.png';

fprintf('1. NIHT, 2. CGIHT, 3. ASD, 4. ScaledASD (suggested:))\n')
alg = input('Select an algorithm number from 1 to 4: ');

fprintf('1. boat, 2. barb, 3. lena\n')
img = input('Select an image number from 1 to 3: ');

fprintf('\n')

mask = 'random';
if strcmp(mask,'random')
  random_inpaint(alg_list{alg},img_list{img})
else
  cross_inpaint(alg_list{alg},img_list{img})
end