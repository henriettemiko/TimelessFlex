# TimelessFlex - A flexible framework for investigating chromatin state trajectories

by Henriette Miko (henriette.miko@mdc-berlin.de), [Ohler lab](
https://github.com/ohlerlab) at BIMSB/MDC, October 03, 2020


## Guix information

TimelessFlex_manifest.scm is the manifest for the GNU Guix profile that was used for TimelessFlex.

A docker container image was created with the following command "guixr pack -m TimelessFlex_manifest.scm -f docker -C none -S /bin=bin -S /lib=lib -S /share=share --save-provenance" and can be downloaded at https://bimsbstatic.mdc-berlin.de/ohler/henriettemiko/TimelessFlex-docker-pack.tar .


For reproducibility the output of the two guix commands "guix describe" and "guix describe -f channels" are given in guix_describe_out.txt and guix_describe_channels_out.txt respectively.



