FROM ubuntu:18.04

RUN apt-get update  && apt-get upgrade && apt-get install -y build-essential

RUN apt install python3.6 -y
RUN apt install python3-pip -y

RUN pip3 install numpy scipy matplotlib ipython jupyter pandas sympy nose astropy peakutils

RUN apt install python3-pyqt5 -y  
RUN apt install pyqt5-dev-tools -y
RUN apt install qttools5-dev-tools -y

RUN pip3 install pyqt5
RUN apt install gedit kate -y

RUN mkdir -p /root/.config/matplotlib
RUN echo "backend : Agg" > /root/.config/matplotlib/matplotlibrc







