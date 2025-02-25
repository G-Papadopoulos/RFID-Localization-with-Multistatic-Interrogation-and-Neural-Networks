Run "rand_data.m" to produce the dataset (for random tag positions), which will include some complex DoA values. These rows will be removed by running "remove_complex_samples_csv_converter.m". This will also create a .csv file to use in Colab.

Many variables in "rand_data.m" can be altered, to create different channels and different topologies.

"grid_data.m" is similar to "rand_data.m", with the positions of the tag mimicking a non-perfect path, e.g. a line with small fluctuations (used for testing).

"Ambiguity" and "Ambiguity (new system)" contain different setups.

Experimental data available only for the original setup.