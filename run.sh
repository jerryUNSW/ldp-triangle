#!/bin/bash

epsilon=2

data="../graphs/email.edges"

# data="../graphs/test2.edges"

# data="../graphs/small.edges"

round=10

num_threads=1

algo=1

./abcore $epsilon $data $round $num_threads $algo