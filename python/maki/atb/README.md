# AllTheBacteria Interface

#### Count assemblies

```{console}
# Species
maki.atb count "Escherichia coli"

# Genus
maki.atb count "Bacillus" --rank genus
```

#### List samples

```{console}
# Basic
maki.atb list "Campylobacter jejuni"

# High-quality assemblies
maki.atb list "Campylobacter jejuni" --high-quality --assemblies

# Limit output
maki.atb list "Escherichia coli" --limit 100
```

#### Plan downloads

```{console}
maki.atb download-plan "Escherichia coli"

# Strict filtering
maki.atb download-plan "Bacillus" \
    --rank genus \
    --high-quality \
    --assemblies
```