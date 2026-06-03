FROM feiler98/utiltosica

RUN mkdir -p /scratch/tmp/feiler/TOSICA_PrepData_cellXgene
WORKDIR /scratch/tmp/feiler/TOSICA_PrepData_cellXgene

COPY . .

# . . means from our computer into our container
RUN pip install --no-cache-dir -r requirements.txt
CMD ["python3", "/scratch/tmp/feiler/TOSICA_PrepData_cellXgene/atlas_prep_tosica.py"]