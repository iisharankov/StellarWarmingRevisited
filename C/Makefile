CC = gcc
LD = gcc   
OPT = -O3 -fopt-info-vec#-missed
CFLAGS = -Wall -g $(OPT)
# INCLUDES = -I/opt/X11/include
LDFLAGS = -lm -lpthread  #-L/opt/X11/lib -lX11
# LDLIBS =

RM = rm -f
OBJS = main.o const.o equations.o starGenerator.o generateHRD.o
NAME = starCreator

all: $(NAME)

$(NAME): $(OBJS)
	$(LD) $(OPT) $(OBJS) $(LDFLAGS)  -o $(NAME)

main.o: main.c const.h equations.h starGenerator.h generateHRD.h Makefile
	$(CC) $(CFLAGS) -c main.c

equations.o: equations.h equations.c
	$(CC) $(CFLAGS) -c equations.c

generateHRD.o: generateHRD.h generateHRD.c
	$(CC) $(CFLAGS) -c generateHRD.c


starGenerator.o: starGenerator.h starGenerator.c
	$(CC) $(CFLAGS) -c starGenerator.c

const.o: const.c const.h
	$(CC) $(CLAGS) -c const.c

clean:
	$(RM) $(NAME) $(OBJS) result.gal
