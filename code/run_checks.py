"""Single fail-closed entry point; execute original assertions even under -O.

Scripts run in temporary directories, leaving published evidence untouched.
Finite exact checks are not a proof of the unrestricted theorem.
"""
import ast
import hashlib
import json
from pathlib import Path
import sys
import tempfile

class PreserveAssertions(ast.NodeTransformer):
    def visit_Assert(self, node):
        return ast.copy_location(ast.If(test=ast.UnaryOp(op=ast.Not(), operand=node.test), body=[ast.Raise(exc=ast.Call(func=ast.Name(id='AssertionError', ctx=ast.Load()), args=[node.msg or ast.Constant(value='Check failed')], keywords=[]), cause=None)], orelse=[]), node)

def execute(path, replacement=None):
    text = path.read_text()
    if replacement:
        old, new = replacement
        if old not in text:
            raise RuntimeError('Mutation anchor missing')
        text = text.replace(old, new, 1)
    tree = ast.fix_missing_locations(PreserveAssertions().visit(ast.parse(text)))
    with tempfile.TemporaryDirectory(prefix='knese-check-') as tmp:
        target = Path(tmp) / path.name
        target.write_text(text)
        exec(compile(tree, str(target), 'exec', optimize=0), {'__name__':'__main__', '__file__':str(target)})

root = Path(__file__).resolve().parent
paths = [root/'exact-operator-checks.py', root/'independent-exact-checks.py', root/'fresh-targeted-checks.py']
for path in paths:
    execute(path)
    print(json.dumps({'checker':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'passed':True}))
mutation = ('weight=F(1,n*choose(n,l))', 'weight=F(1,(n+1)*choose(n,l))')
try:
    execute(paths[0], mutation)
except AssertionError:
    print(json.dumps({'negative_control':'incorrect inverse-binomial weight','rejected':True,'optimized':bool(sys.flags.optimize)}))
else:
    raise RuntimeError('Corrupt weights escaped detection')
print(json.dumps({'all_passed':True,'optimized':bool(sys.flags.optimize),'scope':'Bounded exact consistency checks; no formal or external verification'}))
