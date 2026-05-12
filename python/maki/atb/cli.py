import typer
from pathlib import Path
from .config import DEFAULT_DB
from .db import connect, ensure_db
from .schema import inspect_schema
from .columns import resolve_columns
from .queries import count as count_impl, list_samples
from .download import plan_download as plan_download_impl

app = typer.Typer()


def load(db_path: Path):
    db_path = ensure_db(db_path)
    conn = connect(db_path)

    schema = inspect_schema(conn)
    mapping = resolve_columns(schema)

    return conn, schema, mapping


@app.command()
def init(db_path: Path = DEFAULT_DB):
    """
    Download and prepare metadata database.
    """
    ensure_db(db_path)
    typer.echo(f"Database ready at {db_path}")


@app.command()
def inspect(db_path: Path = DEFAULT_DB):
    """
    Inspect schema.
    """
    conn = connect(db_path)
    schema = inspect_schema(conn)

    import json
    typer.echo(json.dumps(schema, indent=2))


@app.command()
def count(
    taxon: str,
    rank: str = "species",
    db_path: Path = DEFAULT_DB,
):
    conn, schema, mapping = load(db_path)

    result = count_impl(conn, mapping, taxon, rank)
    typer.echo(result)


@app.command()
def list(
    taxon: str,
    rank: str = "species",
    high_quality: bool = False,
    assemblies: bool = False,
    limit: int = 20,
    db_path: Path = DEFAULT_DB,
):
    conn, schema, mapping = load(db_path)

    rows = list_samples(conn, mapping, taxon, rank, high_quality, assemblies, limit)

    for r in rows:
        typer.echo(r)


@app.command()
def plan_download(
    taxon: str,
    rank: str = "species",
    high_quality: bool = True,
    assemblies: bool = True,
    db_path: Path = DEFAULT_DB,
):
    conn, schema, mapping = load(db_path)

    results = plan_download_impl(
        conn, schema, mapping, taxon, rank, high_quality, assemblies
    )

    for archive, count in results:
        typer.echo(f"{archive} ({count} samples)")


def main():
  app()
  

if __name__ == "__main__":
    main()
