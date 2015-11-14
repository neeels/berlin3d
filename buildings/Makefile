CPUs = 4

src-names:
	mkdir src-names
	ln -s ../src-citygml/Charlottenburg-Wilmersdorf/ src-names/cburg_wilmersdf
	ln -s ../src-citygml/Friedrichshain-Kreuzberg/ src-names/fhain_xberg
	ln -s ../src-citygml/Lichtenberg/ src-names/lichtenbg
	ln -s ../src-citygml/Marzahn-Hellersdorf/ src-names/marzahn_hellersdf
	ln -s ../src-citygml/Mitte/ src-names/mitte
	ln -s ../src-citygml/Neukoelln/ src-names/neukoelln
	ln -s ../src-citygml/Pankow/ src-names/pankow
	ln -s ../src-citygml/Reinickendorf/ src-names/reinickendf
	ln -s ../src-citygml/Spandau/ src-names/spandau
	ln -s ../src-citygml/Steglitz-Zehlendorf/ src-names/steglitz_zehlendf
	ln -s ../src-citygml/Tempelhof-Schoeneberg/ src-names/tempelhf_schoenebg
	ln -s ../src-citygml/Treptow-Koepenick/ src-names/treptow_koepenick

.PHONY: new-gml-sections gml-sections new-cpp-sections

new-cpp-sections:
	rm -rf cpp-sections
	make cpp-sections

cpp-sections:
	./prepare_cpp_sections_makefile.sh
	make -j $(CPUs) -C cpp-sections

new-gml-sections:
	rm -rf gml-sections
	make gml-sections

gml-sections:
	mkdir gml-sections
	make -j $(CPUs) gml-section-jobs

gml-section-jobs: gml0 gml1 gml2 gml3 gml4 gml5 gml6 gml7 gml8 gml9 gml10 gml11


gml0:
	./section_citygml.py src-citygml/Charlottenburg-Wilmersdorf

gml1:
	./section_citygml.py src-citygml/Friedrichshain-Kreuzberg

gml2:
	./section_citygml.py src-citygml/Lichtenberg

gml3:
	./section_citygml.py src-citygml/Marzahn-Hellersdorf

gml4:
	./section_citygml.py src-citygml/Mitte

gml5:
	./section_citygml.py src-citygml/Neukoelln

gml6:
	./section_citygml.py src-citygml/Pankow

gml7:
	./section_citygml.py src-citygml/Reinickendorf

gml8:
	./section_citygml.py src-citygml/Spandau

gml9:
	./section_citygml.py src-citygml/Steglitz-Zehlendorf

gml10:
	./section_citygml.py src-citygml/Tempelhof-Schoeneberg

gml11:
	./section_citygml.py src-citygml/Treptow-Koepenick

