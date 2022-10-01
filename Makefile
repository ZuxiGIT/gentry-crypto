server: ./server/server.cpp ./src/*.cpp ./include/*.hpp
	g++ ./server/server.cpp ./src/*.cpp ./include/*.hpp -o server

client: ./client/client.cpp ./src/*.cpp ./include/*.hpp
	g++ ./client/client.cpp ./src/*.cpp ./include/*.hpp -o client