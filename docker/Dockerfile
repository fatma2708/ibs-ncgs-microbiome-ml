# Étape ML : reproductible, sans les gros outils de séquençage
FROM python:3.12-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY src/ src/
COPY workflow/ workflow/
COPY tests/ tests/
COPY data/ data/
ENTRYPOINT ["python", "src/run_ml.py"]
