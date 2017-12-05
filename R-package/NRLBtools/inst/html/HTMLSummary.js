/* JavaScript for REDUCE Suite v2: Xiang-Jun Lu <xl2134@columbia.edu> */

/* get all tables in this document -- global variable */
tables = document.getElementsByTagName('table');

function alternateRowColors() {
    if (!document.getElementById || !tables) {
        return;
    }

    for (var tableNumber = 0; tableNumber < tables.length; tableNumber++) {
        var table = tables[tableNumber];
        var tbodies = table.getElementsByTagName('tbody');

        for (var bodyNumber = 0; bodyNumber < tbodies.length; bodyNumber++) {
            var tbody = tbodies[bodyNumber];
            rows = tbody.getElementsByTagName('tr');
            var even = false;

            for (var rowNumber = 0; rowNumber < rows.length; rowNumber++) {
                var row = rows[rowNumber];

                row.onmouseover = function() {
                    this.className += " ruled";
                    return false;
                }
                row.onmouseout = function() {
                    this.className = this.className.replace("ruled", "");
                    return false
                }

                if (even) {
                    row.className += " even";
                } else {
                    row.className += " odd";
                }
                even = !even;
            }
        }
    }
}

/* the following are global variables */
var lastTableNumber;
var lastColumnNumber;
var lastOrder;

function sortTableRows() {
    if (!document.getElementById || !tables) {
        return;
    }

        /* iterate through all the tables, to make the tables sortable */
    for (var tableNumber = 0; tableNumber < tables.length; tableNumber++) {
        var table = tables[tableNumber];
        var thead = table.getElementsByTagName('thead');
        if (thead.length < 1) {  /* the table must have a <thead>-element */
            continue;
        }

            /* assuming one table has one <thead> with one <tr> */
        var headerRow = thead[0].getElementsByTagName('tr')[0];

            /* iterate through all the th-elements, to add the onclick
             * event with the 'sortColumn'-function attached */
        var columns = headerRow.getElementsByTagName('th');
        for (var columnNumber = 0; columnNumber < columns.length; columnNumber++) {
            var column = columns[columnNumber];
            var attrValue = 'sortColumn(' + tableNumber + ', ' + columnNumber + ');';
            column.setAttribute('onclick', attrValue);
        }
    }
}

function sortColumn(tableNumber, columnNumber) {
        /* check wheter it should be ordered ascending or descending */
    var order = 'ascending';
    if (lastTableNumber == tableNumber && lastColumnNumber == columnNumber
        && lastOrder == 'ascending')
        order = 'descending';

    var table = tables[tableNumber];
    var tbodies = table.getElementsByTagName('tbody');
    if (tbodies.length < 1) {  /* there must be <tbody>-element in the table */
        return;
    }

    for (var bodyNumber = 0; bodyNumber < tbodies.length; bodyNumber++) {
        var tbody = tbodies[bodyNumber];
        var rows = tbody.getElementsByTagName('tr');

            /* iterate through the rows, and put the data according to
             * which the table must be sorted in an array */
        var rowData = new Array();
        for (var rowNumber = 0; rowNumber < rows.length; rowNumber++) {
            var row = rows[rowNumber];
            var data = row.getElementsByTagName('td')[columnNumber].firstChild.nodeValue;
            rowData[rowData.length] = new Array(rowNumber, data);
        }

            /* choose the relevant sorting method, according to the first row: */
        var rowDataSorted = new Array();
        if (isValidNumber(rowData[0][1])) {  /* - the data is numeric */
            rowDataSorted = rowData.sort(numeric);
        } else {
            rowDataSorted = rowData.sort(normal);  /* default to standard sorting-method */
        }

            /* iterate through the rows, that now are in the correct order ... */
        var rowToBeRemoved = new Array();
        var rowToBeInserted = new Array();
        for (var rowNumber = 0; rowNumber < rowDataSorted.length; rowNumber++) {
            var originalRow = rows[rowDataSorted[rowNumber][0]];
                /* create a clone of the orignal row, and register it
                 * to be inserted in its new place ... */
            rowToBeInserted[rowToBeInserted.length] = originalRow.cloneNode(true);
                /* and register the original row to be removed,
                 * because its place isn't correct anymore ... */
            rowToBeRemoved[rowToBeRemoved.length] = originalRow;
        }

            /* reverse order if we sort this column descending */
        if (order == 'descending') {
            rowToBeInserted.reverse();
        }

            /* insert the clones into the table's body */
        var even = false;
        for (var rowNumber = 0; rowNumber < rowToBeInserted.length; rowNumber++) {
            var row = rowToBeInserted[rowNumber];

            row.onmouseover = function() {  /* repeat from alternateRowColors() */
                this.className += " ruled";
                return false;
            }
            row.onmouseout = function() {
                this.className = this.className.replace("ruled", "");
                return false
            }

                /* the following two lines are for the alternate row color */
            var oldClass = row.getAttribute('class');
            oldClass = oldClass.replace('odd', '');
            oldClass = oldClass.replace('even', '');

            if (even) {
                row.setAttribute('class', oldClass + ' even');
            } else {
                row.setAttribute('class', oldClass + ' odd');
            }
            even = !even;

            tbody.appendChild(row);
        }

            /* remove the original, now superfluous, rows */
        for (var rowNumber = 0; rowNumber < rowToBeRemoved.length; rowNumber++) {
            tbody.removeChild(rowToBeRemoved[rowNumber]);
        }
    }

        /* register a class to the header of the column according to
         * which the table has been sorted, so one can style a
         * ascending and descending sorted column differently */
    var thead = table.getElementsByTagName('thead');
    var headerRow = thead[0].getElementsByTagName('tr')[0];

        /* iterate through all the th-elements, to add the 'order' class */
    var columns = headerRow.getElementsByTagName('th');
    for (var columnNumberTwo = 0; columnNumberTwo < columns.length; columnNumberTwo++) {
        var column = columns[columnNumberTwo];
        if (columnNumberTwo == columnNumber) {
            column.setAttribute('class', order);
        } else {
            column.setAttribute('class', '');
        }
    }

        /* remember according to which column from which table we've
         * sorted, and in what order, so when one sorts this column
         * again, we sort it in reverse order */
    lastTableNumber = tableNumber;
    lastColumnNumber = columnNumber;
    lastOrder = order;
}

function numeric(foo, bar) {
    return foo[1] - bar[1];
}

function normal(foo, bar) {
    if (foo[1] < bar[1])
        return -1;
    else if (foo[1] > bar[1])
        return 1;
    else
        return 0;
}

function isValidNumber(inpString) {
    return /^[-+]?\d+(\.\d+)?([eE][-+]\d+)?$/.test(inpString);
}

window.onload = function () {
    sortTableRows();
    alternateRowColors();
}
