class QueryBuilder:
    def __init__(self, core_table):
        self.core = core_table
        self.joins = []
        self.conditions = []
        self.params = []

    def join(self, table, on_left, on_right):
        self.joins.append(
            f"JOIN {table.name} ON {on_left} = {on_right}"
        )
        return self

    def where(self, condition, *params):
        self.conditions.append(condition)
        self.params.extend(params)
        return self

    def build(self, select_clause):
        query = f"SELECT {select_clause} FROM {self.core.name}"

        if self.joins:
            query += " " + " ".join(self.joins)

        if self.conditions:
            query += " WHERE " + " AND ".join(self.conditions)

        return query, self.params
      