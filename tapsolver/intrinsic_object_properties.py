from define_adspecies import define_adspecies
from define_gas import define_gas
from display_gasses import display_gasses
from display_surface import display_surface
from experimental_data import experimental_data
from mechanism import mechanism
from reactor import reactor
from reactor_species import reactor_species
from read_old_input import read_old_input
from TAPobject import TAPobject

#testGen1 = load_example('inl_reactor_1.obj')
testGen1 = readCSV_reactor('./input_file_2.csv')
testGen2 = readCSV_mechanism('./input_file_2.csv')
testGen3 = readCSV_ic('./input_file_2.csv')