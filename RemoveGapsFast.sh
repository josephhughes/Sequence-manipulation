#!/bin/bash

sed '/^>/! s/\-//g' $1 > $2