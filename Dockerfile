FROM python:2

# set a directory for the app
WORKDIR /usr/src/temptrend

# install dependencies
RUN pip install --no-cache-dir numpy matplotlib

# copy all the files to the container
COPY . .

# run the command
# ENTRYPOINT ["./temptrend.py"]
CMD ["/bin/bash"]