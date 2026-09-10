CC = gcc
CFLAGS = -O3 -fopenmp -shared -fPIC
TARGET = libcrt_filter_31_omp.so
SRC = crt_filter_31_omp.c
LIBS = -lgmp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LIBS)

clean:
	rm -f $(TARGET)
