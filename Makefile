VERSION=1.09
FILES=corset.cc Cluster.cc Transcript.cc MakeClusters.cc Read.cc MakeClusters.h Transcript.h Read.h StringSet.h Cluster.h corset_fasta_ID_changer.cc
OBJ=corset.o Cluster.o MakeClusters.o Transcript.o Read.o
EX=corset corset_fasta_ID_changer

all: $(EX) Makefile

corset: $(OBJ)
	g++ -I/group/bioi1/nadiad/software/samtools-0.1.18/ -O3 -ffast-math -I./  -I/usr/local/installed/glibc/2.14/include -I/usr/local/installed/zlib/1.2.8/include -I/usr/local/installed/curl/7.49.1/include -I/usr/local/installed/gcc/4.9.4/include -I/usr/local/installed/gcc/4.9.4/include/c++/4.9.4 -I/usr/local/installed/java/1.8.0_66/include -I/usr/local/installed/java/1.8.0_66/include/linux -std=c++0x -DUNORDEREDMAP $(OBJ) -o $@ -L/group/bioi1/nadiad/software/samtools-0.1.18/ -L/usr/local/installed/glibc/2.14/lib -L/usr/local/installed/zlib/1.2.8/lib -L/usr/local/installed/curl/7.49.1/lib -L/usr/local/installed/gcc/4.9.4/lib64 -L/usr/local/installed/gcc/4.9.4/lib -lbz2 -lcurl -llzma -lbam -lz 

%.o: %.cc
	g++ -I/group/bioi1/nadiad/software/samtools-0.1.18/ -O3 -ffast-math -I./  -I/usr/local/installed/glibc/2.14/include -I/usr/local/installed/zlib/1.2.8/include -I/usr/local/installed/curl/7.49.1/include -I/usr/local/installed/gcc/4.9.4/include -I/usr/local/installed/gcc/4.9.4/include/c++/4.9.4 -I/usr/local/installed/java/1.8.0_66/include -I/usr/local/installed/java/1.8.0_66/include/linux -std=c++0x -DUNORDEREDMAP -c $< -o $@ -DVERSION=$(VERSION)

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
