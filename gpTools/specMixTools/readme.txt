- getGPSMSpec.m returns the spectrum of a spectral mixture kernel and
its derivatives wrt parameters
- sampleGPSM returns a function drawn from a gp with zero mean and a SM
kernel of given hypers
- getObjSMGP returns neg log lik and its derivative, evaluated in
frequency domain using fft
- trainSMGP_freq hypeparameter learning by minimising neg log lik,
in freq domain
- denoiseSMGP_freq denoises given signal, returns posterior mean and
marginal variances, given the hyperparameters
- initSMParams initialises hypers for training


tests/ test code
