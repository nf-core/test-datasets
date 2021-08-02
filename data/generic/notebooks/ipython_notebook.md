```python
import matplotlib.pyplot as plt
```

```python tags=["parameters"]
input_dir = None
artifact_dir = None
cpus = 1
n_iter = 10
input_filename = None
```

## Running from nextflow
The parameters are:

```python
print(f"input_dir = {input_dir}")
print(f"artifact_dir = {artifact_dir}")
print(f"cpus = {cpus}")
print(f"n_iter = {n_iter}")
print(f"input_filename = {input_filename}")
```

## Read an input file

```python
if input_dir is not None:
    with open(f"{input_dir}/{input_filename}") as f:
        print(f.read())
```

## Print Hello world `n` times

```python
for i in range(n_iter):
    print(f"Hello World {i+1}!")
```

## Include a plot

```python
plt.plot(range(10))
```

## Write an output file

```python
if artifact_dir is not None:
    with open(f"{artifact_dir}/artifact.txt", 'w') as f:
        f.write("Hello World!\n")
```

