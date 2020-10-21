"""
Contains the BootstrapForm class.
"""
import flask_wtf
import jinja2
import wtforms as wtf








class BootstrapForm(flask_wtf.FlaskForm):
    """
    Detailed description.
    """


    def fieldsToHTML(
        self
        ):
        """
        Detailed description.
        """
        ret = []
        for field in self:
            ret.append("<div class='form-group'>")
            if field.type == "BooleanField":
                ret += self.__booleanField(field)
            ret.append("<label for='"+field.id+"'>"+field.label.text+"</label>")
            if field.type == "SelectField":
                ret += self.__selectField(field)
            elif field.type == "DecimalField":
                ret += self.__decimalField(field)
            if field.errors:
                for error in field.errors:
                    ret.append("<div class='invalid-feedback'>"+error+"</div>")
            else:
                ret.append("<small class='form-text text-muted'>"+field.description+"</small>")
            ret.append("</div>")
        return jinja2.Markup("\n".join(ret) + "\n")


    def __booleanField(
        self
        ,field
        ):
        """
        Detailed description.

        Parameters
        ----------
        field : object
                Detailed description.
        """
        line = (
            "<input class='form-check-input"
            + self.__validity(field)
            + "' name='"
            + field.name
            + "' id='"
            + field.id
            + "' type='checkbox'"
        )
        if field.data:
            line += " value='y'>"
        else:
            line += " value='n'>"
        return [line]


    def __decimalField(
        self
        ,field
        ):
        """
        Detailed description.

        Parameters
        ----------
        field : object
                Detailed description.
        """
        line = (
            "<input class='form-control"
            + self.__validity(field)
            + "' name='"
            + field.name
            + "' id='"
            + field.id
            + "' type='number'"
        )
        if field.raw_data:
            line += " value='"+field.raw_data[0]+"'>"
        else:
            line += " value='"+str(field.data)+"'>"
        return [line]


    def __selectField(
        self
        ,field
        ):
        """
        Detailed description.

        Parameters
        ----------
        field : object
                Detailed description.
        """
        ret = []
        line = (
            "<select class='form-control"
            + self.__validity(field)
            + "' name='"
            + field.name
            + "' id='"
            + field.id
            + "'>"
        )
        ret.append(line)
        for (value,display) in field.choices:
            line = "<option value='"+value+"'"
            if value == field.data:
                line += " selected>"
            else:
                line += ">"
            line += display+"</option>"
            ret.append(line)
        ret.append("</select>")
        return ret


    def __validity(
        self
        ,field
        ):
        """
        Detailed description.

        Parameters
        ----------
        field : object
                Detailed description.
        """
        if self.errors:
            if field.errors:
                return " is-invalid"
            else:
                return " is-valid"
        else:
            return ""
