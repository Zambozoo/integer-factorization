# integer-factorization
A simple GNU BigNum approach to Integer Factorization. Requires GNU. Limited since it needs to generate all primes <= _n_ where n is the largest factorizable integer.

Build: 
  - gcc -o factor factor.c -lgmp -lm -O3
  - ./factor -build [TBL_DIR="table/"] 
  
Run:
  - ./factor -factor [TRIES=30 [TBL_DIR="table/"]]] 

Basic idea:
- ∃ set of numbers **S**<sub>n</sub> such that for any distinct primes _a_, _b_ ≤ √(_n_), ∃ integer _m_ ∈ **S**<sub>n</sub> such that _a_ | _m_ /\ ~(_b_ | _m_)
- Let |**S**<sub>n</sub>| = _l_, then (_l_ choose (_l_ / 2)) ≥ **pi**(√(_n_)), where **pi**(_x_) = (# of primes <= _x_)
- The min_cardinality(**S**<sub>n</sub>) == _k_, where (_k_ choose _k_ / 2) ≥ **pi**(_n_) /\ ((_k_ - 2) choose (_k_ - 2) / 2) < **pi**(_n_)

- ∃ **S**<sub>n</sub> for all _n_ such that for every _m_ in **S**<sub>n</sub>, all prime factors of _m_ have multiplicity 1 /\ **num_prime_factors**(_m_) == **pi**(√(_n_)) / 2

- Let **T**<sub>n</sub> be the set of all **S**<sub>n</sub> with minimal cardinality, then |**T**<sub>n</sub>| = (|**S**<sub>n</sub>| choose **pi**(_n_))

TLDR; After generating a set of integers **S**<sub>n</sub> where its entries are products of combinations of primes ≤ √(_n_), you are guaranteed to factor any number _m_ ≤ _n_ if GCD(_m_, _i_) > 1 and i ∈ **S**<sub>n</sub>.
