# A Robust RELAX-Based Algorithm for Enhanced Electric Network Frequency Estimation

We would like to acknowledge (Hua et al.) for their work [Robust ENF Estimation Based on Harmonic Enhancement and Maximum Weight Clique](https://ieeexplore.ieee.org/abstract/document/9494518) which inspired us to develop our RELAX-based estimation schemes. We would like to thank them for sharing their codes in [ENF_Enhancement_Estimation](https://github.com/ghua-ac/ENF-WHU-Dataset/tree/master/ENF_Enhancement_Estimation) folder. Moreover, we would like to thank (Li et al.) for inspiring us to employ the relaxation algorithm [Efficient mixed-spectrum estimation with applications to target feature extraction](https://ieeexplore.ieee.org/document/540585) for ENF estimation.

# Our Paper

This repo contains the codes of the proposed RELAX-based ENF estimation described in the [**A Robust RELAX-Based Algorithm for Enhanced Electric Network Frequency Estimation**]() paper. The paper is published in the *** Proceedings of the 13th Conference on Artificial Intelligence*** (SETN 2024).

# Dataset

For the evaluation of our method, we used the H1 and H1_ref folders of the [ENF-WHU-Dataset](https://github.com/ghuawhu/ENF-WHU-Dataset/tree/master/ENF-WHU-Dataset).

# MATLAB Codes

1. **my_relax_test.m** (Code for estimating the ENF for all schemes)
2. **my_plots.m** (Code for generating the figures presented in the paper)
3. **paired_ttest.m** (Code for the hypothesis testing leveraging paired t-tests)

# Citation

```

@inproceedings{korgialas2024relax,
  title={A Robust RELAX-Based Algorithm for Enhanced Electric Network Frequency Estimation},
  author={Korgialas, Christos and Kotropoulos, Constantine},
  booktitle={Proceedings of the 13th Conference on Artificial Intelligence (SETN)},
  pages={1--6},
  year={2024},
  publisher={ACM}
}
```

# Acknowledgements

This research was supported by the Hellenic Foundation for Research and Innovation (H.F.R.I.) under the ``2nd Call for H.F.R.I. Research Projects to support Faculty Members & Researchers" (Project Number: 3888).

# Authors

Feel free to send us a message for any issue.

***Christos Korgialas (ckorgial@csd.auth.gr) and Constantine Kotropoulos (costas@csd.auth.gr)***
