"""
Symbolic processor class for ANDES models.
"""

import logging
import os
import numpy as np
from collections import OrderedDict, defaultdict

from sympy import Symbol, Matrix
from sympy import sympify, lambdify, latex, SympifyError
from sympy import SparseMatrix

from andes.shared import dilled_vars
from andes.utils.paths import get_pycode_path

logger = logging.getLogger(__name__)


class SymProcessor:
    """
    A helper class for symbolic processing and code generation.

    Parameters
    ----------
    parent : Model
        The `Model` instance to document

    Attributes
    ----------
    xy : sympy.Matrix
        variables pretty print in the order of State, ExtState, Algeb, ExtAlgeb
    f : sympy.Matrix
        differential equations pretty print
    g : sympy.Matrix
        algebraic equations pretty print
    df : sympy.SparseMatrix
        df /d (xy) pretty print
    dg : sympy.SparseMatrix
        dg /d (xy) pretty print
    inputs_dict : OrderedDict
        All possible symbols in equations, including variables, parameters, discrete flags, and
        config flags. It has the same variables as what ``get_inputs()`` returns.
    vars_dict : OrderedDict
        variable-only symbols, which are useful when getting the Jacobian matrices.

    """

    def __init__(self, parent):

        self.parent = parent
        # symbols that are input to lambda functions
        # including parameters, variables, services, configs, and scalars (dae_t, sys_f, sys_mva)
        self.inputs_dict = OrderedDict()
        self.lambdify_func = [dict(), 'numpy']

        self.vars_dict = OrderedDict()
        self.vars_int_dict = OrderedDict()   # internal variables only
        self.vars_list = list()

        self.f_list, self.g_list = list(), list()  # symbolic equations in lists
        self.f_matrix, self.g_matrix, self.s_matrix = list(), list(), list()  # equations in matrices

        # pretty print of variables
        self.xy = list()  # variables in the order of states, algebs
        self.f, self.g, self.s = list(), list(), list()
        self.df, self.dg = None, None

        # get references to the parent attributes
        self.calls = parent.calls
        self.cache = parent.cache
        self.config = parent.config
        self.class_name = parent.class_name
        self.tex_names = OrderedDict()

    def generate_symbols(self):
        """
        Generate symbols for symbolic equation generations.

        This function should run before other generate equations.

        Attributes
        ----------
        inputs_dict : OrderedDict
            name-symbol pair of all parameters, variables and configs

        vars_dict : OrderedDict
            name-symbol pair of all variables, in the order of (states_and_ext + algebs_and_ext)

        """

        logger.debug(f'- Generating symbols for {self.class_name}')

        # clear symbols storage
        self.f_list, self.g_list = list(), list()
        self.f_matrix, self.g_matrix = Matrix([]), Matrix([])

        # process tex_names defined in model
        # -----------------------------------------------------------
        for key in self.parent.tex_names.keys():
            self.tex_names[key] = Symbol(self.parent.tex_names[key])
        for instance in self.parent.discrete.values():
            for name, tex_name in zip(instance.get_names(), instance.get_tex_names()):
                self.tex_names[name] = tex_name
        # -----------------------------------------------------------

        for var in self.cache.all_params_names:
            self.inputs_dict[var] = Symbol(var)

        for var in self.cache.all_vars_names:
            tmp = Symbol(var)
            self.vars_dict[var] = tmp
            self.inputs_dict[var] = tmp
            if var in self.cache.vars_int:
                self.vars_int_dict[var] = tmp

        # store tex names defined in `self.config`
        for key in self.config.as_dict():
            tmp = Symbol(key)
            self.inputs_dict[key] = tmp
            if key in self.config.tex_names:
                self.tex_names[tmp] = Symbol(self.config.tex_names[key])

        # store tex names for pretty printing replacement later
        for var in self.inputs_dict:
            if var in self.parent.__dict__ and self.parent.__dict__[var].tex_name is not None:
                self.tex_names[Symbol(var)] = Symbol(self.parent.__dict__[var].tex_name)

        self.inputs_dict['dae_t'] = Symbol('dae_t')
        self.inputs_dict['sys_f'] = Symbol('sys_f')
        self.inputs_dict['sys_mva'] = Symbol('sys_mva')

        self.lambdify_func[0]['Indicator'] = lambda x: x
        self.lambdify_func[0]['imag'] = np.imag
        self.lambdify_func[0]['real'] = np.real
        self.lambdify_func[0]['im'] = np.imag
        self.lambdify_func[0]['re'] = np.real

        self.vars_list = list(self.vars_dict.values())  # useful for ``.jacobian()``

    def _check_expr_symbols(self, expr):
        """
        Check if expression contains unknown symbols.
        """
        fs = expr.free_symbols
        for item in fs:
            if item not in self.inputs_dict.values():
                raise ValueError(f'{self.class_name} expression "{expr}" contains unknown symbol "{item}"')

        return fs

    def generate_equations(self):
        logger.debug(f'- Generating equations for {self.class_name}')

        self.f_list, self.g_list = list(), list()

        self.calls.f = None
        self.calls.g = None
        self.calls.f_args = list()
        self.calls.g_args = list()

        vars_list = [self.cache.states_and_ext, self.cache.algebs_and_ext]
        expr_list = [self.f_list, self.g_list]

        eqn_names = ['f', 'g']
        eqn_args = [self.calls.f_args, self.calls.g_args]

        for vlist, elist, ename, eargs in zip(vars_list, expr_list, eqn_names, eqn_args):
            sym_args = list()
            for name, instance in vlist.items():
                if instance.e_str is None:
                    elist.append(0)
                else:
                    try:
                        expr = sympify(instance.e_str, locals=self.inputs_dict)
                    except (SympifyError, TypeError) as e:
                        logger.error('Error parsing equation "%s "for %s.%s',
                                     instance.e_str, instance.owner.class_name, name)
                        raise e

                    free_syms = self._check_expr_symbols(expr)
                    for s in free_syms:
                        if s not in sym_args:
                            sym_args.append(s)
                            eargs.append(str(s))

                    elist.append(expr)
            if len(elist) == 0 or not any(elist):  # `any`, not `all`
                self.calls.__dict__[ename] = None
            else:
                self.calls.__dict__[ename] = lambdify(sym_args, tuple(elist),
                                                      modules=self.lambdify_func)

        # convert to SymPy matrices
        self.f_matrix = Matrix(self.f_list)
        self.g_matrix = Matrix(self.g_list)

    def generate_services(self):
        """
        Generate calls for services, including ``ConstService`` and
        ``VarService`` among others.

        Sequence is preserved due to possible dependency
        """

        s_args = OrderedDict()
        s_syms = OrderedDict()
        s_calls = OrderedDict()

        for name, instance in self.parent.services.items():
            v_str = '0' if instance.v_str is None else instance.v_str
            try:
                expr = sympify(v_str, locals=self.inputs_dict)
            except (SympifyError, TypeError) as e:
                logger.error(f'Error parsing equation for {instance.owner.class_name}.{name}')
                raise e
            self._check_expr_symbols(expr)
            s_syms[name] = expr
            s_args[name] = [str(i) for i in expr.free_symbols]
            s_calls[name] = lambdify(s_args[name], s_syms[name], modules=self.lambdify_func)

        self.s_matrix = Matrix(list(s_syms.values()))
        self.calls.s = s_calls
        self.calls.s_args = s_args

    def generate_jacobians(self):
        """
        Generate Jacobians and store to corresponding triplets.

        The internal indices of equations and variables are stored, alongside the lambda functions.

        For example, dg/dy is a sparse matrix whose elements are ``(row, col, val)``, where ``row`` and ``col``
        are the internal indices, and ``val`` is the numerical lambda function. They will be stored to

            row -> self.calls._igy
            col -> self.calls._jgy
            val -> self.calls._vgy

        """
        logger.debug(f'- Generating Jacobians for {self.class_name}')

        # clear storage
        self.df_syms, self.dg_syms = Matrix([]), Matrix([])
        self.calls.clear_ijv()

        # NOTE: SymPy does not allow getting the derivative of an empty array
        if len(self.g_matrix) > 0:
            self.dg_syms = self.g_matrix.jacobian(self.vars_list)

        if len(self.f_matrix) > 0:
            self.df_syms = self.f_matrix.jacobian(self.vars_list)

        self.df_sparse = SparseMatrix(self.df_syms)
        self.dg_sparse = SparseMatrix(self.dg_syms)

        vars_syms_list = list(self.vars_dict)
        algebs_and_ext_list = list(self.cache.algebs_and_ext)
        states_and_ext_list = list(self.cache.states_and_ext)

        fg_sparse = [self.df_sparse, self.dg_sparse]
        j_args = defaultdict(list)   # argument list for each jacobian call
        j_calls = defaultdict(list)  # jacobian functions (one for each type)

        for idx, eq_sparse in enumerate(fg_sparse):
            for item in eq_sparse.row_list():
                e_idx, v_idx, e_symbolic = item
                if idx == 0:
                    eq_name = states_and_ext_list[e_idx]
                else:
                    eq_name = algebs_and_ext_list[e_idx]

                var_name = vars_syms_list[v_idx]
                eqn = self.cache.all_vars[eq_name]    # `BaseVar` that corr. to the equation
                var = self.cache.all_vars[var_name]   # `BaseVar` that corr. to the variable
                jname = f'{eqn.e_code}{var.v_code}'

                # jac calls with all arguments and stored individually
                self.calls.append_ijv(jname, e_idx, v_idx, 0)

                # collect unique arguments for jac calls
                free_syms = self._check_expr_symbols(e_symbolic)
                for fs in free_syms:
                    if fs not in j_args[jname]:
                        j_args[jname].append(fs)
                j_calls[jname].append(e_symbolic)

        for jname in j_calls:
            self.calls.j_args[jname] = [str(i) for i in j_args[jname]]
            self.calls.j[jname] = lambdify(j_args[jname], tuple(j_calls[jname]), modules=self.lambdify_func)

        self.calls.j_names = list(j_calls.keys())

        # The for-loop below is intended to add an epsilon small value to the diagonal of `gy`.
        # The user should take care of the algebraic equations by using `diag_eps` in `Algeb` definition

        for var in self.parent.cache.vars_int.values():
            if var.diag_eps == 0.0:
                continue
            elif var.diag_eps is True:
                eps = self.parent.system.config.diag_eps
            else:
                eps = var.diag_eps

            if var.e_code == 'g':
                eq_list = algebs_and_ext_list
            else:
                eq_list = states_and_ext_list

            e_idx = eq_list.index(var.name)
            v_idx = vars_syms_list.index(var.name)

            self.calls.append_ijv(f'{var.e_code}{var.v_code}c', e_idx, v_idx, eps)

    def generate_pretty_print(self):
        """
        Generate pretty print variables and equations.
        """
        logger.debug(f"- Generating pretty prints for {self.class_name}")

        # equation symbols for pretty printing
        self.f, self.g = Matrix([]), Matrix([])

        self.xy = Matrix(list(self.vars_dict.values())).subs(self.tex_names)

        # get pretty printing equations by substituting symbols
        self.f = self.f_matrix.subs(self.tex_names)
        self.g = self.g_matrix.subs(self.tex_names)
        self.s = self.s_matrix.subs(self.tex_names)

        # store latex strings
        nx = len(self.f)
        ny = len(self.g)
        self.calls.x_latex = [latex(item) for item in self.xy[:nx]]
        self.calls.y_latex = [latex(item) for item in self.xy[nx:nx + ny]]

        self.calls.f_latex = [latex(item) for item in self.f]
        self.calls.g_latex = [latex(item) for item in self.g]
        self.calls.s_latex = [latex(item) for item in self.s]

        self.df = self.df_sparse.subs(self.tex_names)
        self.dg = self.dg_sparse.subs(self.tex_names)

        # store init latex strings
        init_latex = OrderedDict()
        for name, instance in self.cache.all_vars.items():
            if instance.v_str is None and instance.v_iter is None:
                init_latex[name] = ''
            else:
                if instance.v_str is not None:
                    init_latex[name] = latex(self.v_str_syms[name].subs(self.tex_names))
                if instance.v_iter is not None:
                    init_latex[name] = latex(self.v_iter_syms[name].subs(self.tex_names))

        self.calls.init_latex = init_latex

    def generate_pycode(self, pycode_path):
        """
        Create output source code file for generated code.

        Generated code are stored at ``~/.andes/pycode``.
        """

        import pprint
        pycode_path = get_pycode_path(pycode_path, mkdir=True)

        file_path = os.path.join(pycode_path, f'{self.class_name}.py')
        header = \
            """from collections import OrderedDict  # NOQA

from numpy import nan, pi, sin, cos, tan, sqrt, exp, select  # NOQA
from numpy import greater_equal, less_equal, greater, less   # NOQA
from numpy import logical_and, logical_or, logical_not  # NOQA
from numpy import array, real, imag, conj, angle, arctan, radians  # NOQA
from numpy import log  # NOQA


"""

        with open(file_path, 'w') as f:
            f.write(header)
            f.write(self._rename_func(self.calls.f, 'f_update'))
            f.write(self._rename_func(self.calls.g, 'g_update'))

            for name in self.calls.j:
                f.write(self._rename_func(self.calls.j[name], f'{name}_update'))

            # initialization: assignments
            for name in self.calls.ia:
                f.write(self._rename_func(self.calls.ia[name], f'{name}_ia'))
            for name in self.calls.ii:
                f.write(self._rename_func(self.calls.ii[name], f'{name}_ii'))
            for name in self.calls.ij:
                f.write(self._rename_func(self.calls.ij[name], f'{name}_ij'))

            # services
            for name in self.calls.s:
                f.write(self._rename_func(self.calls.s[name], f'{name}_svc'))

            # variables
            for name in dilled_vars:
                f.write(f'\n{name} = ' + pprint.pformat(self.calls.__dict__[name]))

    def _rename_func(self, func, func_name):
        """
        Rename the function name and return source code.

        This function does not check for name conflicts.
        Install `yapf` for optional code reformatting (takes extra processing time).
        """
        import inspect

        if func is None:
            return f"# empty {func_name}\n"

        src = inspect.getsource(func)
        src = src.replace("def _lambdifygenerated(", f"def {func_name}(")
        # remove `Indicator`
        src = src.replace("Indicator", "")

        if self.parent.system.config.yapf_pycode:
            try:
                from yapf.yapflib.yapf_api import FormatCode
                src = FormatCode(src, style_config='pep8')[0]  # drop the encoding `None`
            except ImportError:
                logger.warning("`yapf` not installed. Skipped code reformatting.")

        src += '\n'
        return src

    def generate_dependency(self):
        """
        Generate dependency list and initialization order.
        """
        self.v_str_syms = OrderedDict()
        self.v_iter_syms = OrderedDict()
        deps = dict()

        # convert to symbols
        for name, instance in self.cache.all_vars.items():
            if instance.v_str is not None:
                sympified = sympify(instance.v_str, locals=self.inputs_dict)
                self._check_expr_symbols(sympified)
                self.v_str_syms[name] = sympified
            else:
                # default initial values to zero
                sympified = sympify('0.0', locals=self.inputs_dict)
                self.v_str_syms[name] = sympified

            if instance.v_iter is not None:
                sympified = sympify(instance.v_iter, locals=self.inputs_dict)
                self._check_expr_symbols(sympified)
                self.v_iter_syms[name] = sympified

        # store deps for explicit and iterative initializers
        for name, expr in self.v_str_syms.items():
            _store_deps(name, expr, self.vars_dict, deps)

        for name, expr in self.v_iter_syms.items():
            _store_deps(name, expr, self.vars_dict, deps)

        # resolve dependency
        self.init_seq = resolve_deps(deps)
        self.calls.init_seq = self.init_seq

    def check_v_iter(self):
        """
        Helper function to check if `v_iter` is defined for variables
        with circular dependencies.
        """

        for item in self.init_seq:
            if not isinstance(item, list):
                continue

            for vi in item:
                if self.cache.all_vars[vi].v_iter is None:
                    logger.error("%s: v_iter not defined for %s" % (self.class_name, vi))

    def generate_init(self):
        """
        Generate initialization equations.
        """

        self.generate_dependency()
        self.check_v_iter()
        self.generate_init_eqn()
        self.lambdify_init()

    def generate_init_eqn(self):
        """
        Generate initialization equations.

        The RHS of assignment equations ``v = v_str(x, y)`` or RHS of iterative
        initialization equations in the form of ``0 = v_iter(x, y)`` will be
        stored to ``self.init_list``.

        A list of flags will be stored to ``self.init_flag`` with 0 for
        assignments and 1 for iterative.

        For iteratively initialized variables that require assigned initial
        values (to improve convergence), the initial value can be provided
        through `self.v_str`.
        """

        self.init_asn = OrderedDict()       # assignment-type initialization
        self.init_itn = OrderedDict()       # iterative initialization
        self.init_itn_vars = OrderedDict()  # variables corr. to iterative vars
        self.init_jac = OrderedDict()

        for item in self.init_seq:
            if isinstance(item, str):
                instance = self.parent.__dict__[item]
                if instance.v_str is not None:
                    self.init_asn[item] = self.v_str_syms[item]
                if instance.v_iter is not None:
                    self.init_itn[item] = self.v_iter_syms[item]

            elif isinstance(item, list):
                name_concat = '_'.join(item)
                eqn_set = Matrix([self.v_iter_syms[name] for name in item])
                self.init_itn[name_concat] = eqn_set
                self.init_itn_vars[name_concat] = item
                for vv in item:
                    instance = self.parent.__dict__[vv]
                    if instance.v_str is not None:
                        self.init_asn[vv] = self.v_str_syms[vv]

        for name, expr in self.init_itn.items():
            vars_iter = OrderedDict()
            for item in self.init_itn_vars[name]:
                vars_iter[item] = self.vars_dict[item]

            self.init_jac[name] = expr.jacobian(list(vars_iter.values()))

    def lambdify_init(self):
        """
        Convert equations and Jacobians to lambda functions.
        """

        init_a = OrderedDict()
        init_i = OrderedDict()
        init_j = OrderedDict()
        ia_args = OrderedDict()  # arguments for assignment init.
        ii_args = OrderedDict()  # arguments for iterative init.
        ij_args = OrderedDict()

        for name, expr in self.init_asn.items():
            self._check_expr_symbols(expr)
            ia_args[name] = [str(i) for i in expr.free_symbols]
            init_a[name] = lambdify(ia_args[name], expr, modules=self.lambdify_func)

        for name, expr in self.init_itn.items():
            self._check_expr_symbols(expr)
            ii_args[name] = [str(i) for i in expr.free_symbols]
            init_i[name] = lambdify(ii_args[name], expr, modules=self.lambdify_func)

            jexpr = self.init_jac[name]
            ij_args[name] = [str(i) for i in jexpr.free_symbols]
            init_j[name] = lambdify(ij_args[name], jexpr, modules=self.lambdify_func)

        self.calls.ia = init_a
        self.calls.ii = init_i
        self.calls.ij = init_j
        self.calls.ia_args = ia_args
        self.calls.ii_args = ii_args
        self.calls.ij_args = ij_args


def _store_deps(name, sympified, vars_int_dict, deps):
    """
    Helper function to store dependencies to a dict.

    Used by ``resolve``.
    """

    deps[name] = []
    for fs in sympified.free_symbols:
        if fs not in vars_int_dict.values():
            continue
        if fs not in deps[name]:
            deps[name].append(str(fs))


def resolve_deps(graph):
    """
    Resolve dependency for a dict-based graph using recursion.
    """

    seq = list()      # sequence after resolving dependency
    visited = list()  # book keeper of the visited nodes

    cflat = list()    # flattened node in circles
    circles = list()  # circles as lists

    def sub_resolve(name, deps, path):
        for item in deps:
            if (item == name) or (item in visited):
                continue

            # undeclared leaf node
            if (item not in graph):
                if item not in seq:
                    seq.append(item)
                continue

            # circular dependency
            if item in path:
                idx1 = path.index(item)
                idx2 = path.index(name)
                mn = min(idx1, idx2)
                mx = max(idx1, idx2)
                cflat.extend(path[mn:mx+1])
                circles.append(path[mn:mx+1])
                continue

            path.append(item)
            sub_resolve(item, graph[item], path)

        # when all dependent nodes are visited
        if name not in visited:
            # if the current node is not in any circle
            if name not in cflat:
                seq.append(name)

            else:
                for cc in circles:
                    if (name in cc) and (cc not in seq):
                        seq.append(cc)

            visited.append(name)

    for name, deps in graph.items():
        path = list()
        sub_resolve(name, deps, path)

    return seq
