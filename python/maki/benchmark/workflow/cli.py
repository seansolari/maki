
from typing import Annotated

import typer


Workflow = Annotated[
    str,
    typer.Option("-w", "--workflow",
                 help="Workflow to run")
]
