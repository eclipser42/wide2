TARGET_DIR = build/WildlifeDensity\ gfortran.build/Development/WildlifeDensity.build/Objects-normal/ppc

$(TARGET_DIR)/mnps2.o: mnps2.f
	gfortran -g -fno-underscoring -c -o "$@" $<
