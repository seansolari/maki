
from typing import Annotated

import typer


Workflow = Annotated[
    str,
    typer.Option("-w", "--workflow",
                 help="Workflow to run")
]


def parse_workflow_args(ctx: typer.Context):
    extra_options = {}
    it = iter(ctx.args)
    for item in it:
        if item.startswith("-"):
            key = item.lstrip("-").replace("-", "_")
            extra_options[key] = next(it, True)
    return extra_options
