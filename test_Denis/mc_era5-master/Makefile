# !!!IMPORTANT!!!
# make should be run in the appropriate python (conda) environment!
# CONDA_ACTIVATE=conda activate clim
# CONDA_DEACTIVATE=conda deactivate
EXEC_NB=python execute_notebook.py

CODE_DIR=code
FIG_DIR=figures
DATA_DIR=data

FIGURES=\
    $(FIG_DIR)/ascat_era5_interim_accacia_case_vort_wspd.pdf \
    $(FIG_DIR)/vrf__vort_thresh__tfreq__bs2000_100.pdf \
    $(FIG_DIR)/characteristic_histograms.pdf \
    $(FIG_DIR)/density_maps__track_genesis_lysis.pdf

DATA_IN=\
    $(DATA_DIR)/tracks/stars/PolarLow_tracks_North_2002_2011

ASCAT=$(DATA_DIR)/"ascat"/"ascat_20130326_123601_metopa_33384_eps_o_125_2101_ovw.l2.nc"

AVHRR=$(DATA_DIR)/"avhrr"/"201303261220.ch4.r8.nps.tif"
AVHRR_CASSINI=$(DATA_DIR)/"avhrr"/"201303261220.ch4.r8.cass.tif"

$(AVHRR): $(AVHRR_CASSINI)

$(FIG_DIR)/ascat_era5_interim_accacia_case_vort_wspd.pdf: $(CODE_DIR)/ACCACIA-Case-Example.ipynb $(AVHRR) $(ASCAT)
$(FIG_DIR)/vrf__vort_thresh__tfreq__bs2000_100.pdf: $(CODE_DIR)/Verification.ipynb
$(FIG_DIR)/characteristic_histograms.pdf: $(CODE_DIR)/Characteristics.ipynb
$(FIG_DIR)/density_maps__track_genesis_lysis.pdf: $(CODE_DIR)/Density-Maps.ipynb


all: $(FIGURES)

$(FIGURES):
	@echo "making figure: $@"
	$(EXEC_NB) $<

# Reproject the AVHRR GeoTIFF image to a North Pole Stereo projection
$(AVHRR):
	gdalwarp $< $@ -t_srs "+proj=stere +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs"

.PHONY: clean help
clean:
	@echo "Cleaning..."
	rm -f $(FIG_DIR)/*

help:
	@echo ""
	@echo "Usage:"
	@echo "    make all: run Jupyter Notebooks to create all figures"
	@echo "    make clean: Delete files in figures/ folder"
	@echo "    make help: Print this message and exit"
	@echo ""
