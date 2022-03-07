## Supported Polynomial Commitment Schemes

The library supports four polynomial commitment schemes.

### Inner-product-argument PC

A polynomial commitment scheme based on the hardness of the discrete logarithm problem in prime-order groups. 

The construction is described in the following paper.

[pcd-acc]: https://ia.cr/2020/499

[Proof-Carrying Data from Accumulation Schemes][pcd-acc]     
Benedikt Bünz, Alessandro Chiesa, Pratyush Mishra, Nicholas Spooner     
TCC 2020

### Marlin variant of the Kate-Zaverucha-Goldberg PC

[kzg10]: http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf
[marlin]: https://ia.cr/2019/1047

Polynomial commitment based on the Kate-Zaverucha-Goldberg construction, with degree enforcement, batching, and (optional) hiding property taken from Marlin.

The construction is described in the following paper.

[Marlin: Preprocessing zkSNARKs with Universal and Updatable SRS][marlin]     
Alessandro Chiesa, Yuncong Hu, Mary Maller, Pratyush Mishra, Noah Vesely, Nicholas Ward  
EUROCRYPT 2020

[Polynomial Commitments][kzg10]     
Aniket Kate, Gregory M. Zaverucha, Ian Goldberg     
ASIACRYPT 2010

### Sonic/AuroraLight variant of the Kate-Zaverucha-Goldberg PC

Polynomial commitment based on the Kate-Zaverucha-Goldberg construction, with degree enforcement and batching taken from Sonic (more precisely, their counterparts in AuroraLight that avoid negative G1 powers). The (optional) hiding property of the commitment scheme follows the approach described in Marlin.

The construction is described in the following papers.

[sonic]: https://ia.cr/2019/099
[aurora-light]: https://ia.cr/2019/601

[AuroraLight: Improved Prover Efficiency and SRS Size in a Sonic-Like System][aurora-light]     
Ariel Gabizon     
ePrint, 2019

[Sonic: Zero-Knowledge SNARKs from Linear-Size Universal and Updateable Structured Reference Strings][sonic]     
Mary Maller, Sean Bowe, Markulf Kohlweiss, Sarah Meiklejohn     
CCS 2019

[Marlin: Preprocessing zkSNARKs with Universal and Updatable SRS][marlin]     
Alessandro Chiesa, Yuncong Hu, Mary Maller, Pratyush Mishra, Noah Vesely, Nicholas Ward  
EUROCRYPT 2020

[Polynomial Commitments][kzg10]     
Aniket Kate, Gregory M. Zaverucha, Ian Goldberg     
ASIACRYPT 2010

### Marlin variant of the Papamanthou-Shi-Tamassia multivariate PC

Multivariate polynomial commitment based on the construction in the Papamanthou-Shi-Tamassia construction with batching and (optional) hiding property inspired by the univariate scheme in Marlin.

The construction is described in the following paper.

[pst]: https://ia.cr/2011/587

[Signatures of Correct Computation][pst]    
Charalampos Papamanthou, Elaine Shi, Roberto Tamassia   
TCC 2013

[Marlin: Preprocessing zkSNARKs with Universal and Updatable SRS][marlin]     
Alessandro Chiesa, Yuncong Hu, Mary Maller, Pratyush Mishra, Noah Vesely, Nicholas Ward  
EUROCRYPT 2020

## Comparison

### Comparison of `MarlinKZG10` and `SonicKZG10`


#### High-level:
They handle degree bounds differently. 

MarlinPC uses shift powers only in G1 and requires two commitments to enforce degree bounds.

SonicPC uses shift powers in G1 and G2 and requires only one commitment to enforce degree bounds.

#### Setup:

SonicPC additionally computes some G2 elements for shift powers: `(1/\beta)^i H`. This results in a longer verifying key, as shift powers in SonicPC are in G2, while shift powers in Marlin are in G1, and are shared with the "non-shift" powers.

#### Commit:

When there is no degree bound, both are the same.

When there is a degree bound, MarlinPC is more expensive: it needs an additional commitment to commit to the shifted poynomial. 

#### Open: 

When there is no degree bound, both are the same.

When there is a degree bound, MarlinPC is slightly more expensive: it requires more scalar field computations.

#### Check:

MarlinPC simply adjusts the commitment of the shifted polynomial, so the overhead is small. It checks a pairing equation with two pairing operations.

SonicPC is more expensive, as it checks a pairing equation of three pairing operations. It can be reduced into two if there is no degree bound.
