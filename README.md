# RTSMS
Randomized Tucker with single-mode sketching

This repository contains an implementation of RTSMS which is designed to approximate the Tucker decomposition of a given tensor of order d. In the landscape of randomized tensor decompositions, it is common practice to use sketching from multiple sides simultaneously. A distinct feature of RTSMS is that it takes the simpler approach of sketching from just one side at a time. The methodology draws inspiration from the generalized Nystrom and the sequentially truncated HOSVD while incorporating leverage scores, regularization, and iterative refinement to ensure both efficiency and stability in the computation of factor matrices.

Users have the flexibility to input either their desired multilinear rank or a tolerance on the residual. In the latter case, RTSMS automatically determines the appropriate multilinear rank.

Requirements
=============================

Tensorlab is needed to run the codes in this repository. Tensorlab is available at https://www.tensorlab.net. 


Reference
=============================
For a detailed explanation of RTSMS and its applications, please see the following manuscript:

B. Hashemi and Y. Nakatsukasa, RTSMS: randomized Tucker with single-mode sketching, submitted, 2023, https://arxiv.org/abs/2311.14873.
