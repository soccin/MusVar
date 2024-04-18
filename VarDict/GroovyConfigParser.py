from shlex import shlex
from ast import literal_eval

TRANSLATION = {
    "true": True,
    "false": False,
    "null": None,
}


class ParseException(Exception):
    def __init__(self, token, line):
        self.token = token
        self.line = line

    def __str__(self):
        return "ParseException at line %d: invalid token %s" % (self.line, self.token)


class GroovyConfigParser:
    def __init__(self, source):
        if isinstance(source, str):
            self.source = open(source)
            self.should_close_source = True
        else:
            self.source = source
            self.should_close_source = False

    def __del__(self):
        if self.should_close_source and not self.source.closed:
            self.source.close()

    def parse(self):
        lex = shlex(self.source)
        lex.wordchars = lex.wordchars + ".-"
        lex.commenters = "//"
        state = 1
        context = []
        result = dict()
        while True:
            token = lex.get_token()
            if not token:
                return result
            if state == 1:
                if token == "}":
                    if len(context):
                        context.pop()
                    else:
                        raise ParseException(token, lex.lineno)
                else:
                    name = token
                    state = 2
            elif state == 2:
                if token == "=":
                    state = 3
                elif token == "{":
                    context.append(name)
                    state = 1
                else:
                    raise ParseException(token, lex.lineno)
            elif state == 3:
                try:
                    value = TRANSLATION[token]
                except KeyError:
                    value = literal_eval(token)
                key = ".".join(context + [name]).split(".")
                current = result
                for i in range(len(key) - 1):
                    if key[i] not in current:
                        current[key[i]] = dict()
                    current = current[key[i]]
                current[key[-1]] = value
                state = 1

    def write_as_json_file(self, json_file):
        import json
        with open(json_file, 'w') as file:
            json.dump(self.parse(), file, indent=4)

