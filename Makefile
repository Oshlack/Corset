VERSION=1.07
FILES=corset.cc Cluster.cc Transcript.cc MakeClusters.cc Read.cc MakeClusters.h Transcript.h Read.h StringSet.h Cluster.h corset_fasta_ID_changer.cc
OBJ=corset.o Cluster.o MakeClusters.o Transcript.o Read.o
EX=corset corset_fasta_ID_changer

all: $(EX) Makefile

corset: $(OBJ)
	g++ -I/group/bioi1/nadiad/software/samtools-0.1.18 -O3 -ffast-math -I./   -std=c++0x -DUNORDEREDMAP $(OBJ) -o $@ -L/group/bioi1/nadiad/software/samtools-0.1.18  -lpthread -lbam -lz 

%.o: %.cc
	g++ -I/group/bioi1/nadiad/software/samtools-0.1.18 -O3 -ffast-math -I./   -std=c++0x -DUNORDEREDMAP -c $< -o $@ -DVERSION=$(VERSION)

corset_fasta_ID_changer: corset_fasta_ID_changer.cc
	g++ $< -o $@

clean:                                                     
	-rm *~ $(OBJ)

Clean: clean
	-rm $(EX)

install: all
	test -d /usr/local/bin/ || mkdir /usr/local/bin/
	cp $(EX) /usr/local/bin/

tar_ball:
	-mkdir corset-$(VERSION)
	cp $(FILES) corset-$(VERSION)
	cp Makefile.in README configure COPYING LICENSE corset-$(VERSION)
	tar -cvf corset-$(VERSION).tar corset-$(VERSION)/
	gzip corset-$(VERSION).tar

tar_ball_binary:
	-mkdir corset-$(VERSION)-linux64
	cp $(EX) corset-$(VERSION)-linux64
	cp COPYING LICENSE corset-$(VERSION)-linux64
	tar -cvf corset-$(VERSION)-linux64.tar corset-$(VERSION)-linux64/
	gzip corset-$(VERSION)-linux64.tar
