#!/bin/bash

g++ -std=c++11 -pthread src/main.cpp -lboost_program_options -lboost_filesystem -lboost_system
