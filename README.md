Authors: Yu ICHIDA, Rika YAMADA, Shiori KATO, Yuki KAMAYA, Minami KOSUGE, Mamoru AIZAWA, Takashi Okuda SAKAMOTO, Shigetoshi YAZAKI

Title: A simple mathematical model for evaluation of non-fragmentation property of injectable calcium-phosphate cement

Journal: Scientific Reports

parameter settings in the paper

select_phi=1: only on the outer boundary: phi_int+phi_init+RND

RK, L=1, N=200, D0=0.01, tmax=1, T0=0.8, eps=0.1, dt=0.1/N^2, h=2L/N, c: IPc=(0, 0.02, 0.030, phi_bar=0.5

How to proceed: ./doALL.sh

Parameters can be changed in the following files: doALL.sh, evolution/head_paste.h

Shell script files (bash): do.sh, doALL.sh

Source files in evolution/paste.c, head_paste.c, head_paste.h, paste-mkdir.c

Source files in tools/matutil.c, matutil.h

Save data files in data/

Save movie files in movie/

Save backup files in backup/
