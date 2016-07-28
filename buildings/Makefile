CPUs = 4

all: cpp-sections gml-sections

.PHONY: new-gml-sections new-cpp-sections

new-cpp-sections:
	rm -rf cpp-sections
	make cpp-sections

cpp-sections: gml-sections
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


test_cpp_sections: all
	g++ -std=c++11 -Iinclude -Icpp-sections -o test_cpp_sections test_cpp_sections.cpp cpp-sections/*.o

clean-yes-im-sure:
	rm -rf gml-sections cpp-sections
