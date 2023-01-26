
PLUGIN_NAME =  rtxi_module_kernik

HEADERS = rtxi_module_kernik.hpp

SOURCES = rtxi_module_kernik.cpp kernik19.cpp ./include/RealTimeMath.cpp ./include/PowFast.cpp

LIBS =

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile

clean:
	rm -f $(OBJECTS)
	rm -f moc_*
	rm -f *.o
	rm -f $(PLUGIN_NAME).la
	rm -f $(PLUGIN_NAME).o
	rm -rf .libs
	rm -rf include/.libs
	rm -f include/*.o
	rm -f include/moc_*
