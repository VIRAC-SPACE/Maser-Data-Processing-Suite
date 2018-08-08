FROM python:3

RUN pip3 install numpy scipy matplotlib ipython jupyter pandas sympy nose astropy peakutils

# WORKDIR /DataProcessingForMaserObservation
ADD ./code/ /
ADD ./config/ /config
ARG SOURCE

CMD [ "python3", "-u", "./main.py", "cepa"]


