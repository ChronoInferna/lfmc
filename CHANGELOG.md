## 0.2.0 (2026-03-03)

### Feat

- finished refactor and added control variates
- payoffs and estimator
- estimator and payoff outline
- manager config and small tweaks
- antithetic variates + refactor
- **Payoff**: new payoff concept
- **Simulator.hpp**: new abstraction of simulator which uses its own thread
- more constraints on concepts for clarity and compile-time sizing of number of testing threads
- no variance reduction strategy and small updates to polymorphic constructors
- **ThreadPool.hpp**: threadpool for managing threads specific to monte carlo simulations
- **manager.hpp**: create framework for stochastic processes and numerical schemes

### Fix

- fix vr strategy concept
- small work on simulator stuff

### Refactor

- sort cpp files into folders and remove archive
- renamed everything, finished pipeline, and converted each step to return vectors of results
- another huge refactor breaking things into different engines
- simplify simulator interface
- refactor some of xander's code
- remove threadpool
- i don't even know lol
- turns out we don't need to the changes i made to the concepts a little bit ago
- **timing.hpp**: rename from timer
- **NumericalScheme.hpp**: adjust concept for schemes for more freedom
- compile-time strategies for process and scheme and runtime strategy for variance reduction
- flesh out the architecture a bit more
- **Timer**: extremely small return change

## 0.1.0 (2025-12-18)

### Feat

- **ScopedTimer**: raii timer
- **timer.cpp**: created timer

### Refactor

- **timer.cpp**: change method names
- remove codeql (#2)
