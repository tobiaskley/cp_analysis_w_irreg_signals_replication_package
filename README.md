This folder contains the files for simulations and empirical example as used for the revision to AoS in July 2024.

To replicate the tables and figures in the manuscript.
The simulations were performed on the HPC of GWDG in Goettingen, which uses SLURM.
Proceed with the following steps to replicate the results from scratch.

(1A) simulate means and variances of (Z'_t):
* Run /1A_sim_mean_var/simulation_mean.sh
  -> this generates files res_X_bar_[1-300].Rdata
* Run /1A_sim_mean_var/readat.R
  -> this generates Z_mean_var.Rdata
  !! the generated files need to be in folder /res_files

(1B) simulate discretized Brownian bridges
* Run /1B_sim_Q_bridge/simulation_Q_bridge.sh
  -> this generates files res_Q_bridge_[1-50].Rdata
* Run /1B_sim_Q_bridge/readat.R
  -> this generates Q_Bridge.Rdata

!! The files Z_mean_var.Rdata and Q_Bridge.Rdata
   are available in the replication package.

(2A) Real Data Analysis
* Run /2A_real_data/covid.R
  -> creates Figure 1 --- /outPlots/Baidu_index_data.pdf
  -> creates Figure 2 --- /outPlots/coug_I_seq.pdf
                          /outPlots/fev_I_seq.pdf
  -> creates Figure 6 --- /outPlots/covid_comp.pdf
  -> numbers mentioned in the text (Section 5) are returned at the console

(2Ba) Simulation results test
* Run /2Ba_sim_test/sim_test.sh
  -> this generates files result_alt_[1-2000].Rdata
* Run /2Ba_sim_test/readat.R
  -> contents for Tables 2 are returned at the console
  -> creates Figure 4 --- /outPlots/pow_plot.pdf

(2Bb) Simulation results test with long-run variance estimate
* Run /2Bb_sim_test_est_lrv/testing_sim.sh
  -> this generates files result_alt_[1-2000].Rdata
* Run /2Bb_sim_test_est_lrv/readat.R
  -> contents for Tables 3 are returned at the console
  -> creates Figure 7 --- /outPlots/pow_plot_hat_sigma.pdf

(2C) Trend illustration and simulation results locating method
* Run /2C_sim_locating/simulation_alternative.sh
  -> this generates files result_alt_[1-120].Rdata
* Run /2C_sim_locating/readat.R
  -> creates Figure 3 --- /outPlots/illustration_trend.pdf
  -> creates Figure 5 --- /outPlots/bar_plot_five_methods.pdf
  -> creates Figure 8 --- /outPlots/box_plot_five_methods.pdf

