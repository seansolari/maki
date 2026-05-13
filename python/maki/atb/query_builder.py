from copy import deepcopy


class QueryBuilder:
    def __init__(self, core_table):
        self.core = core_table
        self.joins = []
        self.conditions = []
        self.params = []
        self.group_cols = []
        self.order_clause = None
        self.limit = None

    def join(self, table, on_left, on_right):
        self.joins.append(
            f"JOIN {table.name} ON {on_left} = {on_right}"
        )
        return self

    def where(self, condition, *params):
        self.conditions.append(condition)
        self.params.extend(params)
        return self

    def group_by(self, *cols):
        self.group_cols.extend(cols)
        return self

    def order_by(self, clause):
        self.order_clause = clause
        return self
      
    def set_limit(self, n: int):
        self.limit = n

    def build(self, select_clause):
        query = f"SELECT {select_clause} FROM {self.core.name}"
        params = deepcopy(self.params)

        if self.joins:
            query += " " + " ".join(self.joins)

        if self.conditions:
            query += " WHERE " + " AND ".join(self.conditions)

        if self.group_cols:
            query += " GROUP BY " + ", ".join(self.group_cols)

        if self.order_clause:
            query += f" ORDER BY {self.order_clause}"
            
        if self.limit:
            query += " LIMIT ?"
            params.append(self.limit)

        return query, params
      